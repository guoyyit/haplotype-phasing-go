package phase

import (
	"fmt"
	"math"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
	"github.com/guoyyit/haplotype-phasing-go/internal/stats"
)

// adjudicatedEdge is the deterministic decision for one adjacent het pair.
type adjudicatedEdge struct {
	fromHet   int
	connect   bool
	flip      bool // H1 of next het is the opposite bit of previous
	support   int  // votes agreeing with the majority orientation
	conflict  int  // votes disagreeing
	strength  float64
	errorRate float64
	reason    string
}

// rng64 is a deterministic SplitMix64 stream seeded per run; it is only
// consulted for exact evidence ties, never otherwise.
type rng64 struct{ state uint64 }

func (r *rng64) uint64() uint64 {
	r.state += 0x9e3779b97f4a7c15
	z := r.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}

func (r *rng64) bit() int { return int(r.uint64() & 1) }

// Run phases every sample in the dataset deterministically.
func Run(ds *model.Dataset, cfg config.Config) *model.Result {
	res := &model.Result{Samples: map[string]*model.SampleResult{}}
	res.Order = append(res.Order, ds.Samples...)
	res.Stats = model.RunStats{
		Samples:   len(ds.Samples),
		Variants:  len(ds.Variants),
		TrioCount: len(ds.Trios),
	}

	trioByChild := map[string]*model.Trio{}
	for _, trio := range ds.Trios {
		t := trio
		trioByChild[trio.Child] = &t
	}

	for _, sample := range ds.Samples {
		// Each sample derives an independent stream from (seed, sample) so
		// that adding/removing other samples cannot change a sample's ties.
		rng := &rng64{state: mixSeed(cfg.Seed, sample)}
		sr := phaseSample(ds, cfg, sample, trioByChild[sample], rng, &res.Stats)
		res.Samples[sample] = sr
	}

	// Read and trio totals are sample-independent counts over the dataset.
	for _, read := range ds.Reads {
		res.Stats.ReadObservations += len(read.Observations)
	}
	return res
}

// mixSeed derives a per-sample deterministic RNG state from the run seed.
func mixSeed(seed int64, sample string) uint64 {
	h := uint64(seed) ^ 0x9e3779b97f4a7c15
	for i := 0; i < len(sample); i++ {
		h ^= uint64(sample[i])
		h *= 1099511628211
		h ^= h >> 31
	}
	r := &rng64{state: h}
	return r.uint64()
}

// phaseSample phases one sample and accumulates run-level statistics.
func phaseSample(ds *model.Dataset, cfg config.Config, sample string, trio *model.Trio, rng *rng64, runStats *model.RunStats) *model.SampleResult {
	sr := &model.SampleResult{Sample: sample}
	sites := classifySites(ds, sample, cfg.MinSiteDepth)

	// Count missing/partial genotypes and low-coverage sites once per sample.
	for _, gt := range ds.GT[sample] {
		if !gt.Present || gt.Ploidy == 0 {
			runStats.MissingGenotypes++
		}
	}

	// Trio evidence: violations and forced paternal orientation.
	var forced []bool
	var forcedBits []int
	if cfg.EnableTrio && trio != nil {
		forced, forcedBits = trioEvidence(ds, *trio, sites)
	}
	for _, info := range sites {
		if info.hetero {
			sr.HetsTotal++
		}
		if info.violated {
			runStats.TrioViolations++
			if info.hetero {
				sr.ViolationHets++
			}
		}
		if info.lowCov {
			runStats.LowCoverageSites++
		}
	}
	if cfg.EnableTrio && trio != nil {
		runStats.TrioConstraintSites += sr.HetsTotal
	}

	edges := []voteTally{}
	if cfg.EnableReads {
		edges = collectReadEdges(ds, sample, sites)
		runStats.ReadPairs += countReadPairs(edges)
	}
	if cfg.EnableTrio && trio != nil {
		edges = addTrioEdges(edges, forced, forcedBits)
	}

	hetSites := nonViolatedHets(sites)
	decisions := adjudicate(edges, len(hetSites), forced, forcedBits, cfg, rng, runStats)
	blocks := buildBlocks(sample, hetSites, decisions, forced, forcedBits, cfg)

	sr.Sites = buildSitePhases(ds, sample, sites, blocks)
	sr.Blocks = blocks.phaseBlocks
	sr.PhasedHets, sr.UnphasedHets = countPhased(sr.Sites)
	runStats.PhasedHets += sr.PhasedHets
	runStats.UnphasedHets += sr.UnphasedHets
	runStats.HetsTotal += sr.HetsTotal
	runStats.BlocksTotal += len(sr.Blocks)
	return sr
}

// blockLayout maps every non-violated het ordinal to a block id and the
// orientation bit of that het's A0 allele on haplotype H1.
type blockLayout struct {
	phaseBlocks []model.PhaseBlock
	blockOfHet  []int
	bitOfHet    []int
}

// adjudicate applies the documented voting and threshold policy.
func adjudicate(edges []voteTally, hetCount int, forced []bool, forcedBits []int, cfg config.Config, rng *rng64, st *model.RunStats) []adjudicatedEdge {
	byFrom := map[int]int{}
	for i, edge := range edges {
		byFrom[edge.fromHet] = i
	}
	out := make([]adjudicatedEdge, 0, max0(hetCount-1))
	for h := 0; h+1 < hetCount; h++ {
		e := adjudicatedEdge{fromHet: h}
		pos, ok := byFrom[h]
		if !ok {
			e.reason = "NO_EVIDENCE"
			e.strength = 0
			e.errorRate = 0.5
			out = append(out, e)
			continue
		}
		tally := edges[pos]
		// Count a read-vs-trio orientation conflict when both sources cast
		// votes (of any direction) and their directions disagree.
		readVotes := tally.readSame + tally.readOpp
		trioVotes := tally.trioSame + tally.trioOpp
		if readVotes > 0 && trioVotes > 0 {
			agree := (tally.readSame > 0 && tally.trioSame > 0) ||
				(tally.readOpp > 0 && tally.trioOpp > 0)
			oppose := (tally.readSame > 0 && tally.trioOpp > 0) ||
				(tally.readOpp > 0 && tally.trioSame > 0)
			if oppose && !agree {
				st.ReadEdgeConflicts++
			}
		}
		same, opp := tally.totals()
		total := same + opp
		e.support = max0(same)
		e.conflict = opp
		if same == opp {
			// Deterministic, seeded tie-break. The conflict is counted.
			e.flip = rng.bit() == 1
			if e.flip {
				e.support, e.conflict = opp, same
			}
			e.strength = 0.5
			e.errorRate = stats.EdgeErrorRate(e.support, e.conflict)
			st.ReadConflictsAdjudicated++
			st.ReadEdgeConflicts += same
			e.reason = "TIE_SEEDED"
		} else if same > opp {
			e.flip = false
			e.support = same
			e.conflict = opp
			e.strength = float64(same) / float64(total)
			if opp > 0 {
				st.ReadConflictsAdjudicated++
			}
			st.ReadEdgeConflicts += min0(tally.readSame, tally.readOpp)
			if tally.readOpp > 0 && tally.trioSame > 0 {
				st.ReadEdgeConflicts++
			}
			e.errorRate = stats.EdgeErrorRate(same, opp)
			e.reason = "MAJORITY_SAME"
		} else {
			e.flip = true
			e.support = opp
			e.conflict = same
			e.strength = float64(opp) / float64(total)
			if same > 0 {
				st.ReadConflictsAdjudicated++
			}
			st.ReadEdgeConflicts += min0(tally.readSame, tally.readOpp)
			if tally.readSame > 0 && tally.trioOpp > 0 {
				st.ReadEdgeConflicts++
			}
			e.errorRate = stats.EdgeErrorRate(opp, same)
			e.reason = "MAJORITY_OPP"
		}
		e.connect = total >= cfg.MinEdgeSupport &&
			e.strength >= cfg.MinEdgeStrength &&
			e.errorRate <= cfg.BlockMaxErrorRate
		if !e.connect {
			st.BreaksByThreshold++
		} else if e.support > 0 {
			st.ReadEdgeSupports++
		}
		out = append(out, e)
	}
	return out
}

// buildBlocks walks het ordinals, starting a new block on every break.
func buildBlocks(sample string, hetSites []siteInfo, decisions []adjudicatedEdge, forced []bool, forcedBits []int, cfg config.Config) blockLayout {
	layout := blockLayout{}
	n := len(hetSites)
	layout.blockOfHet = make([]int, n)
	layout.bitOfHet = make([]int, n)
	if n == 0 {
		return layout
	}

	// Relative orientation: bitRelative[h] is the edge-implied flip parity.
	relative := make([]int, n)
	for h := 1; h < n; h++ {
		relative[h] = relative[h-1]
		if decisions[h-1].connect && decisions[h-1].flip {
			relative[h] ^= 1
		}
	}

	// Block membership: break on every non-connecting edge.
	blockOf := make([]int, n)
	b := 0
	for h := 1; h < n; h++ {
		if !decisions[h-1].connect {
			b++
		}
		blockOf[h] = b
	}

	// Absolute orientation: if any het in a block has a forced paternal
	// bit, anchor the whole block so that site's H1 carries the paternal
	// allele; otherwise default anchor bit 0 (fully deterministic).
	anchor := map[int]int{}
	if len(forced) == n {
		for h := 0; h < n; h++ {
			if forced[h] {
				blk := blockOf[h]
				want := forcedBits[h]
				// relative[h] xor startBit == want => startBit = ...
				startBit := relative[h] ^ want
				if prev, ok := anchor[blk]; ok && prev != startBit {
					// Conflicting forced bits inside one block: keep the
					// earliest and count via edge conflicts at call sites.
					continue
				}
				anchor[blk] = startBit
			}
		}
	}
	for h := 0; h < n; h++ {
		blk := blockOf[h]
		startBit, ok := anchor[blk]
		if !ok {
			startBit = 0
		}
		layout.bitOfHet[h] = relative[h] ^ startBit
		layout.blockOfHet[h] = blk
	}

	type span struct{ start, end int }
	spans := []span{}
	start := 0
	for h := 1; h <= n; h++ {
		if h == n || blockOf[h] != blockOf[h-1] {
			spans = append(spans, span{start, h - 1})
			start = h
		}
	}
	for blk, spn := range spans {
		start, end := spn.start, spn.end
		var edgeRates []float64
		support, conflict := 0, 0
		for h := start; h < end; h++ {
			d := decisions[h]
			edgeRates = append(edgeRates, d.errorRate)
			support += d.support
			conflict += d.conflict
		}
		errRate := stats.BlockErrorRate(edgeRates)
		confidence := 1 - errRate
		hets := end - start + 1
		status := model.BlockStatusOK
		if hets == 1 {
			status = model.BlockStatusSingleton
			confidence, errRate = 1, 0
		} else if confidence < cfg.MinBlockConfidence {
			status = model.BlockStatusLowConf
		}
		layout.phaseBlocks = append(layout.phaseBlocks, model.PhaseBlock{
			ID:                fmt.Sprintf("BLK-%s-%04d", sample, blk+1),
			StartIdx:          hetSites[start].index,
			EndIdx:            hetSites[end].index,
			SiteCount:         hetSites[end].index - hetSites[start].index + 1,
			HetsCount:         hets,
			Support:           support,
			Conflict:          conflict,
			Confidence:        round4(confidence),
			ExpectedErrorRate: round4(errRate),
			Status:            status,
		})
	}

	return layout
}

// buildSitePhases materializes the per-site rows.
func buildSitePhases(ds *model.Dataset, sample string, sites []siteInfo, layout blockLayout) []model.SitePhase {
	out := make([]model.SitePhase, len(sites))
	ord := map[int]int{}
	h := 0
	for _, info := range sites {
		if info.hetero && !info.violated {
			ord[info.index] = h
			h++
		}
	}
	for i, info := range sites {
		sp := model.SitePhase{
			Index:     i,
			VariantID: ds.Variants[i].ID,
			LowCov:    info.lowCov,
			Violation: info.violated,
			A0:        info.a0,
			A1:        info.a1,
		}
		switch {
		case info.violated:
			sp.Phased = false
			sp.Hetero = info.hetero
			if info.present && info.ploidy == 2 {
				gt := ds.GT[sample][i]
				call := model.AlleleA0
				if gt.A1 != info.a0 && gt.A1 == info.a1 {
					call = model.AlleleA1
				}
				if info.hetero {
					// Record the two observed canonical alleles in the
					// documented unordered form (A0 then A1). The site stays
					// UNPHASED/violation=1, but the original genotype and
					// all violation statistics survive the write/re-read.
					sp.H1, sp.H2 = model.AlleleA0, model.AlleleA1
				} else {
					sp.H1, sp.H2 = call, call
				}
			}
		case info.hetero:
			sp.Phased = true
			sp.Hetero = true
			ordinal := ord[i]
			if layout.bitOfHet[ordinal] == 0 {
				sp.H1, sp.H2 = model.AlleleA0, model.AlleleA1
			} else {
				sp.H1, sp.H2 = model.AlleleA1, model.AlleleA0
			}
		case info.present && info.ploidy == 2:
			gt := ds.GT[sample][i]
			call := model.AlleleA0
			if gt.A1 != info.a0 && gt.A1 == info.a1 {
				call = model.AlleleA1
			}
			sp.H1, sp.H2 = call, call
		case info.present && info.ploidy == 1:
			gt := ds.GT[sample][i]
			if gt.A1 == info.a0 {
				sp.H1 = model.AlleleA0
			} else {
				sp.H1 = model.AlleleA1
			}
		}
		out[i] = sp
	}
	return out
}

// ---- small utilities ----

func countHets(sites []siteInfo) int {
	n := 0
	for _, info := range sites {
		if info.hetero {
			n++
		}
	}
	return n
}

func nonViolatedHets(sites []siteInfo) []siteInfo {
	out := make([]siteInfo, 0)
	for _, info := range sites {
		if info.hetero && !info.violated {
			out = append(out, info)
		}
	}
	return out
}

func countReadPairs(edges []voteTally) int {
	n := 0
	for _, e := range edges {
		n += e.readSame + e.readOpp
	}
	return n
}

func countPhased(sites []model.SitePhase) (phased, unphased int) {
	for _, sp := range sites {
		if !sp.Hetero {
			continue
		}
		if sp.Phased {
			phased++
		} else {
			unphased++
		}
	}
	return phased, unphased
}

func max0(x int) int {
	if x < 0 {
		return 0
	}
	return x
}

func min0(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func round4(v float64) float64 {
	return math.Round(v*10000) / 10000
}

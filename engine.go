package phasing

import (
	"fmt"
	"math"
	"sort"
)

// Run performs the full deterministic phasing pipeline.
func Run(in *Input) (*Result, error) {
	if err := validateInput(in); err != nil {
		return nil, err
	}
	ctx, err := buildContext(in)
	if err != nil {
		return nil, err
	}
	phaseSamples(ctx)
	if in.Trio != nil {
		applyTrio(ctx)
	}
	assembleResult(ctx)
	return ctx.result, nil
}

// Rephase re-reads a prior result and re-emits the identical phasing.
func Rephase(in *Input, prior *Result) (*Result, error) {
	if err := validateInput(in); err != nil {
		return nil, err
	}
	if prior == nil {
		return nil, fmt.Errorf("rephase requires a prior result")
	}
	if prior.SchemaVersion != SchemaVersion {
		return nil, fmt.Errorf("rephase: unsupported schema_version %q", prior.SchemaVersion)
	}
	if err := validatePrior(in, prior); err != nil {
		return nil, err
	}
	return reconstructFromPrior(in, prior), nil
}

// genotypeClass categorizes one diploid call.
type genotypeClass int

const (
	gtMissing genotypeClass = iota // both alleles unknown
	gtPartial                      // exactly one allele known
	gtHomRef
	gtHomAlt
	gtHet
)

func classify(c GenotypeCall) genotypeClass {
	if !c.Present {
		return gtMissing
	}
	switch {
	case c.A1 < 0 && c.A2 < 0:
		return gtMissing
	case c.A1 < 0 || c.A2 < 0:
		return gtPartial
	case c.A1 == 0 && c.A2 == 0:
		return gtHomRef
	case c.A1 == 1 && c.A2 == 1:
		return gtHomAlt
	default:
		return gtHet
	}
}

func validateInput(in *Input) error {
	if in == nil {
		return fmt.Errorf("nil input")
	}
	if len(in.Variants) == 0 {
		return fmt.Errorf("no variant sites provided")
	}
	for i := range in.Samples {
		if len(in.Samples[i].GT) != len(in.Variants) {
			return fmt.Errorf("sample %q has %d genotype calls, expected %d",
				in.Samples[i].ID, len(in.Samples[i].GT), len(in.Variants))
		}
	}
	if in.Config.MinEdgeSupport < 1 {
		return fmt.Errorf("config not resolved: min_edge_support < 1")
	}
	return nil
}

// ---- internal context (populated stage by stage) ----

type sampleCtx struct {
	index    int
	id       string
	gt       []GenotypeCall
	class    []genotypeClass
	hetIdx   []int              // variant indices of heterozygous sites, in canonical order
	coverage []int              // per-variant count of read observations
	edges    map[int]*phaseEdge // key: index into hetIdx of the left het
	dsu      *DSU
	phase    []int8 // per-variant bit for hets: allele placed on h1
	sites    []SitePhase
	blocks   []BlockInfo
}

// phaseEdge summarizes adjudicated votes between two consecutive heterozygous
// sites (in canonical variant order).
type phaseEdge struct {
	hetI, hetJ    int
	votesSame     int
	votesFlip     int
	readIDsSame   []string
	readIDsFlip   []string
	accepted      bool
	parity        int8
	confidence    float64
	conflictReads int
}

type engineCtx struct {
	in              *Input
	samples         map[string]*sampleCtx
	order           []*sampleCtx
	conflict        ConflictReport
	siteAnnotations map[string]map[int]siteAnnotation
	result          *Result
}

func buildContext(in *Input) (*engineCtx, error) {
	ctx := &engineCtx{
		in:              in,
		samples:         make(map[string]*sampleCtx, len(in.Samples)),
		siteAnnotations: make(map[string]map[int]siteAnnotation),
	}
	for i := range in.Samples {
		s := &in.Samples[i]
		sc := &sampleCtx{
			index:    i,
			id:       s.ID,
			gt:       s.GT,
			class:    make([]genotypeClass, len(s.GT)),
			coverage: make([]int, len(s.GT)),
			phase:    make([]int8, len(s.GT)),
		}
		for v, call := range s.GT {
			sc.class[v] = classify(call)
			if sc.class[v] == gtHet {
				sc.hetIdx = append(sc.hetIdx, v)
			}
		}
		sc.edges = make(map[int]*phaseEdge, max0(len(sc.hetIdx)-1))
		ctx.samples[s.ID] = sc
		ctx.order = append(ctx.order, sc)
	}
	ctx.conflict.Entries = []ConflictEntry{}
	ctx.result = &Result{
		SchemaVersion: SchemaVersion,
		Variants:      in.Variants,
		Conflicts:     ctx.conflict,
	}
	return ctx, nil
}

func phaseSamples(ctx *engineCtx) {
	in := ctx.in

	// Read-derived votes are collected per sample. Reads are processed in file
	// order and observation lists are canonically ordered, making collection
	// deterministic.
	for _, rd := range in.Reads {
		sc := ctx.samples[rd.Sample]
		if sc == nil {
			continue // validated during parsing; defensive only
		}

		// Only heterozygous, fully-known sites provide phase information.
		hetObs := make([]SiteObs, 0, len(rd.Obs))
		for _, ob := range rd.Obs {
			sc.coverage[ob.Variant]++
			// An observation disagreeing with a known homozygous genotype is
			// an evidence/genotype conflict: it is counted and cannot vote.
			if ob.Variant >= 0 && ob.Variant < len(sc.class) {
				switch sc.class[ob.Variant] {
				case gtHomRef:
					if ob.Allele != 0 {
						ctx.addConflict(ConflictEntry{
							SampleID:    sc.id,
							ReadID:      rd.ID,
							VariantI:    ob.Variant,
							VariantJ:    ob.Variant,
							VotedParity: int(ob.Allele),
							KeptParity:  0,
							Kind:        "read_genotype",
						})
						continue
					}
				case gtHomAlt:
					if ob.Allele != 1 {
						ctx.addConflict(ConflictEntry{
							SampleID:    sc.id,
							ReadID:      rd.ID,
							VariantI:    ob.Variant,
							VariantJ:    ob.Variant,
							VotedParity: int(ob.Allele),
							KeptParity:  1,
							Kind:        "read_genotype",
						})
						continue
					}
				}
			}
			if ob.Variant >= 0 && ob.Variant < len(sc.class) && sc.class[ob.Variant] == gtHet {
				hetObs = append(hetObs, ob)
			}
		}
		for k := 0; k+1 < len(hetObs); k++ {
			a, b := hetObs[k], hetObs[k+1]
			i, okI := hetPosition(sc, a.Variant)
			j, okJ := hetPosition(sc, b.Variant)
			if !okI || !okJ || j != i+1 {
				continue // no edge between non-consecutive heterozygous sites
			}
			edge := sc.edges[i]
			if edge == nil {
				edge = &phaseEdge{hetI: i, hetJ: j}
				sc.edges[i] = edge
			}
			// A read carries both observed alleles on the same molecule: if
			// they differ the adjacency is "flip", otherwise "same".
			if a.Allele != b.Allele {
				edge.votesFlip++
				edge.readIDsFlip = append(edge.readIDsFlip, rd.ID)
			} else {
				edge.votesSame++
				edge.readIDsSame = append(edge.readIDsSame, rd.ID)
			}
		}
	}

	// Adjudicate every edge deterministically and build parity constraints.
	for _, sc := range ctx.order {
		sc.dsu = NewDSU(len(sc.hetIdx))
		for i := 0; i+1 < len(sc.hetIdx); i++ {
			edge, ok := sc.edges[i]
			if !ok {
				continue
			}
			adjudicateEdge(ctx, sc, edge)
			if edge.accepted {
				sc.dsu.Union(edge.hetI, edge.hetJ, edge.parity)
			}
		}
		assignOrientations(sc)
	}
}

// hetPosition returns the position of variant v in the sample's het list.
func hetPosition(sc *sampleCtx, v int) (int, bool) {
	// hetIdx is sorted by canonical variant order; binary search.
	lo, hi := 0, len(sc.hetIdx)
	for lo < hi {
		mid := (lo + hi) / 2
		switch {
		case sc.hetIdx[mid] < v:
			lo = mid + 1
		case sc.hetIdx[mid] > v:
			hi = mid
		default:
			return mid, true
		}
	}
	return 0, false
}

// adjudicateEdge applies the documented majority/confidence policy, records
// every losing read as a conflict, and marks whether the edge joins a block.
func adjudicateEdge(ctx *engineCtx, sc *sampleCtx, edge *phaseEdge) {
	cfg := ctx.in.Config
	sort.Strings(edge.readIDsSame)
	sort.Strings(edge.readIDsFlip)

	// Majority verdict; ties use the fixed, configurable tie policy. Nothing is
	// ever dropped silently: losing reads are counted below.
	if edge.votesFlip == edge.votesSame {
		edge.parity = resolveTie(ctx, sc, edge)
	} else if edge.votesFlip > edge.votesSame {
		edge.parity = 1
	} else {
		edge.parity = 0
	}

	edge.confidence = edgeConfidence(edge.votesFlip, edge.votesSame, cfg.ReadErrorRate)
	totalVotes := edge.votesFlip + edge.votesSame
	edge.accepted = edge.confidence >= cfg.MinBlockConfidence && totalVotes >= cfg.MinEdgeSupport

	// Record every read voting against the adjudicated parity.
	losing := edge.readIDsSame
	losingParity := 0
	if edge.parity == 0 {
		losing = edge.readIDsFlip
		losingParity = 1
	}
	edge.conflictReads = len(losing)
	for _, id := range losing {
		ctx.addConflict(ConflictEntry{
			SampleID:    sc.sampleID(),
			ReadID:      id,
			VariantI:    sc.hetIdx[edge.hetI],
			VariantJ:    sc.hetIdx[edge.hetJ],
			VotedParity: losingParity,
			KeptParity:  int(edge.parity),
			Kind:        "read_read",
		})
	}
}

// resolveTie resolves an equal-vote edge: seeded PRNG or smallest read id.
func resolveTie(ctx *engineCtx, sc *sampleCtx, edge *phaseEdge) int8 {
	cfg := ctx.in.Config
	if cfg.SeedEnabled {
		// Deterministic per-edge stream: seed mixed with sample/edge position.
		r := newRNG(cfg.RandomSeed + int64(sc.index)*1_000_003 + int64(edge.hetI))
		return int8(r.intn(2))
	}
	// Unseeded: deterministic lexicographic rule over participating reads.
	smallest := ""
	picked := -1
	for _, id := range edge.readIDsSame {
		if picked < 0 || id < smallest {
			smallest, picked = id, 0
		}
	}
	for _, id := range edge.readIDsFlip {
		if picked < 0 || id < smallest {
			smallest, picked = id, 1
		}
	}
	if picked < 0 {
		return 0
	}
	return int8(picked)
}

// assignOrientations sets phase[v] for each het: the allele placed on h1.
// DSU variables are the h1 allele values themselves. Each component is
// anchored on its canonically smallest het-list position (h1 = 0 there),
// which fixes the otherwise arbitrary global flip of every block
// deterministically, independent of DSU root selection.
func assignOrientations(sc *sampleCtx) {
	compIDs := sc.dsu.Components()
	minMember := make(map[int]int)
	rootParity := make(map[int]int8)
	for k := range sc.hetIdx {
		_, parity := sc.dsu.Find(k)
		cid := compIDs[k]
		if cur, seen := minMember[cid]; !seen || k < cur {
			minMember[cid] = k
			rootParity[cid] = parity
		}
	}
	for k, v := range sc.hetIdx {
		_, parity := sc.dsu.Find(k)
		cid := compIDs[k]
		// h1[k] = parity[k] XOR flip, where flip forces h1[min] = 0.
		flip := rootParity[cid]
		sc.phase[v] = parity ^ flip
	}
}

func applyTrio(ctx *engineCtx) {
	in := ctx.in
	trio := in.Trio
	child := ctx.samples[trio.Child]
	father := ctx.samples[trio.Father]
	mother := ctx.samples[trio.Mother]

	annot := make(map[int]siteAnnotation, len(in.Variants))

	// Per-site Mendelian check, independent of read evidence. A violation is
	// recorded for any child site whose alleles cannot be inherited.
	patAlleles := make([]int8, len(in.Variants)) // allele transmitted by father, -1 unknown
	matAlleles := make([]int8, len(in.Variants))
	for v := range in.Variants {
		patAlleles[v], matAlleles[v] = -1, -1
		violation, pat, mat := mendelCheck(child.gt[v], father.gt[v], mother.gt[v])
		a := siteAnnotation{}
		if violation {
			a.mendelViolation = true
		} else if child.class[v] == gtHet {
			// Origin is resolved whenever at least one parent's transmitted
			// allele is pinned; for homozygous child sites both haplotypes
			// carry the same allele anyway.
			if pat >= 0 || mat >= 0 || child.class[v] != gtHet {
				a.originResolved = true
			}
		}
		patAlleles[v], matAlleles[v] = pat, mat
		annot[v] = a
	}
	ctx.siteAnnotations[child.id] = annot

	// Rebuild the child's DSU with pedigree-first priority: trio constraints
	// are unioned first, read-derived edges are added only when consistent.
	// Every contradicted read edge is counted as a read_trio conflict and is
	// never silently dropped.
	child.dsu = NewDSU(len(child.hetIdx))

	// Parent-origin constraints connect child heterozygous sites whose
	// transmitted paternal (or maternal) allele is fixed at both sites.
	// Such constraints can extend blocks beyond read evidence.
	for k := 0; k+1 < len(child.hetIdx); k++ {
		v1, v2 := child.hetIdx[k], child.hetIdx[k+1]
		if annot[v1].mendelViolation || annot[v2].mendelViolation {
			continue
		}
		// A pinned transmitted allele from either parent at both adjacent hets
		// yields a same/flip constraint (the two child haplotypes always carry
		// opposite alleles at a heterozygous site).
		pinned := func(alleles []int8) (int8, bool) {
			if alleles[v1] >= 0 && alleles[v2] >= 0 {
				return alleles[v1] ^ alleles[v2], true
			}
			return 0, false
		}
		if p, ok := pinned(patAlleles); ok {
			child.dsu.Union(k, k+1, p)
		}
		if p, ok := pinned(matAlleles); ok {
			child.dsu.Union(k, k+1, p)
		}
	}

	// Add accepted read edges that do not contradict pedigree constraints.
	for k := 0; k+1 < len(child.hetIdx); k++ {
		edge := child.edges[k]
		if edge == nil || !edge.accepted {
			continue
		}
		v1, v2 := child.hetIdx[k], child.hetIdx[k+1]
		if annot[v1].mendelViolation || annot[v2].mendelViolation {
			continue
		}
		if !child.dsu.Union(k, k+1, edge.parity) {
			_, keptParity := child.dsu.Find(k)
			_, keptParity2 := child.dsu.Find(k + 1)
			ctx.addConflict(ConflictEntry{
				SampleID:    child.id,
				VariantI:    v1,
				VariantJ:    v2,
				VotedParity: int(edge.parity),
				KeptParity:  int(keptParity ^ keptParity2),
				Kind:        "read_trio",
			})
		}
	}

	// Recompute orientations after trio constraints were added.
	assignOrientations(child)

	// Pedigree-resolved singleton hets not connected by read or trio edges
	// still need a definite orientation: h1 is the paternal allele when it is
	// pinned, otherwise the maternal allele (which is the complementary one).
	compIDsAfter := child.dsu.Components()
	compMembers := make(map[int]int)
	for _, cid := range compIDsAfter {
		compMembers[cid]++
	}
	for k, v := range child.hetIdx {
		if compMembers[compIDsAfter[k]] > 1 {
			continue
		}
		a := annot[v]
		if a.mendelViolation {
			continue
		}
		switch {
		case patAlleles[v] >= 0:
			child.phase[v] = patAlleles[v]
		case matAlleles[v] >= 0:
			child.phase[v] = 1 - matAlleles[v]
		}
	}

	// Finalize origin markers with definitive orientations.
	for v, a := range annot {
		if a.mendelViolation || !a.originResolved {
			continue
		}
		if child.class[v] == gtHet {
			// When only the maternal allele is pinned, paternalH1 correctly
			// stays false even though origin is resolved via the mother.
			a.paternalH1 = patAlleles[v] == child.phase[v]
		}
		annot[v] = a
	}
}

// mendelCheck reports whether a child call is inconsistent with its parents,
// and when resolvable returns the specific allele transmitted by each parent.
func mendelCheck(child, father, mother GenotypeCall) (violation bool, paternal, maternal int8) {
	cc := classify(child)
	fc := classify(father)
	mc := classify(mother)
	if cc == gtMissing || cc == gtPartial {
		return false, -1, -1
	}

	// Exact set of alleles a parent can transmit:
	// 0/0->{0}, 1/1->{1}, 0/1->{0,1}, ./{x}->{x}, ./.->unknown {0,1}.
	transmitSet := func(c GenotypeCall, cls genotypeClass) [2]bool {
		var set [2]bool
		switch cls {
		case gtHomRef:
			set[0] = true
		case gtHomAlt:
			set[1] = true
		case gtHet, gtMissing:
			set[0], set[1] = true, true
		case gtPartial:
			known := c.A1
			if known < 0 {
				known = c.A2
			}
			set[known] = true
		}
		return set
	}
	fromFather := transmitSet(father, fc)
	fromMother := transmitSet(mother, mc)

	// Try both parental assignments of the child's two (unordered) alleles.
	try := func(pf, pm int8) bool {
		return fromFather[pf] && fromMother[pm]
	}
	a, b := child.A1, child.A2
	var viable [][2]int8
	if try(a, b) {
		viable = append(viable, [2]int8{a, b})
	}
	if a != b && try(b, a) {
		viable = append(viable, [2]int8{b, a})
	}
	if len(viable) == 0 {
		return true, -1, -1
	}

	// Pin a transmitted allele only when every viable assignment agrees and the
	// transmitting parent could not have passed the other allele.
	pfSet := [2]bool{true, true}
	pmSet := [2]bool{true, true}
	for k, pair := range viable {
		if k == 0 {
			pfSet, pmSet = [2]bool{}, [2]bool{}
		}
		pfSet[pair[0]] = true
		pmSet[pair[1]] = true
	}
	pin := func(set [2]bool, canOther bool) int8 {
		var only int8 = -1
		count := 0
		for allele, ok := range set {
			if ok {
				only, count = int8(allele), count+1
			}
		}
		if count == 1 && !canOther {
			return only
		}
		return -1
	}
	// For a heterozygous child, the parent's transmitted allele is informative
	// only when that parent itself transmits exactly one allele.
	fatherFixed := fc == gtHomRef || fc == gtHomAlt || fc == gtPartial
	motherFixed := mc == gtHomRef || mc == gtHomAlt || mc == gtPartial
	pat := pin(pfSet, !fatherFixed)
	mat := pin(pmSet, !motherFixed)

	// Homozygous child: both parents necessarily transmitted the same allele;
	// pin against fixed parents directly.
	if a == b {
		if fatherFixed {
			pat = a
		}
		if motherFixed {
			mat = a
		}
	}
	return false, pat, mat
}

func assembleResult(ctx *engineCtx) {
	in := ctx.in
	res := ctx.result
	res.Config = ResolvedConfig{Config: in.Config}

	// Conflict counts (entry list is shared; aggregate by deterministic scan).
	for _, e := range ctx.conflict.Entries {
		switch e.Kind {
		case "read_read":
			ctx.conflict.ReadReadConflicts++
		case "read_trio":
			ctx.conflict.ReadTrioConflicts++
		case "read_genotype":
			ctx.conflict.ReadGenotypeConflicts++
		}
	}

	sampleResults := make([]SampleResult, 0, len(ctx.order))
	var stats RunStats
	stats.Samples = len(in.Samples)
	stats.VariantSites = len(in.Variants)
	violatingHets := 0

	for _, sc := range ctx.order {
		sr := SampleResult{ID: sc.id}

		// Components over het-list positions; only adjudicated accepted edges
		// were unioned, so rejected/missing edges already break blocks here.
		compIDs := sc.dsu.Components()
		compSize := make(map[int]int)
		for _, cid := range compIDs {
			compSize[cid]++
		}
		// Block id 0 means "no block". Multi-het components receive ids in the
		// order of their canonically smallest member.
		blockIDByComp := make(map[int]int)
		nextBlock := 1
		for _, cid := range compIDs {
			if compSize[cid] < 2 {
				continue
			}
			if _, assigned := blockIDByComp[cid]; !assigned {
				blockIDByComp[cid] = nextBlock
				nextBlock++
			}
		}

		annotations := ctx.siteAnnotations[sc.id]
		singletonTrioBlock := make(map[int]int) // het-list position -> block id

		sr.Sites = make([]SitePhase, len(in.Variants))
		lowCoverage := make(map[int]bool)
		for v := range sc.class {
			if sc.class[v] != gtMissing && sc.coverage[v] < in.Config.MinReadsPerSite {
				lowCoverage[v] = true
			}
		}

		for v := 0; v < len(in.Variants); v++ {
			site := SitePhase{
				VariantIdx:  v,
				VariantID:   in.Variants[v].ID,
				ReadSupport: sc.coverage[v],
				LowCoverage: lowCoverage[v],
				H1:          -1,
				H2:          -1,
			}
			call := sc.gt[v]
			switch sc.class[v] {
			case gtMissing:
				site.Status = SiteMissing
				stats.MissingSites++
			case gtPartial:
				known := call.A1
				if known < 0 {
					known = call.A2
				}
				site.Status = SiteUnphased
				site.H1, site.H2 = known, -1
				stats.PartialSites++
			case gtHomRef, gtHomAlt:
				allele := int8(0)
				if sc.class[v] == gtHomAlt {
					allele = 1
				}
				site.Status = SitePhased
				site.H1, site.H2 = allele, allele
			case gtHet:
				stats.HeterozygousSites++
				k, _ := hetPosition(sc, v)
				bid := blockIDByComp[compIDs[k]]
				trioAnchored := annotations != nil &&
					annotations[v].originResolved && !annotations[v].mendelViolation
				if bid == 0 && trioAnchored {
					// Single-het pedigree-resolved site: emitted as a
					// singleton trio block (no read evidence required).
					bid = nextBlock
					singletonTrioBlock[k] = bid
					nextBlock++
				}
				if bid == 0 {
					site.Status = SiteUnphased
					site.H1, site.H2 = call.A1, call.A2
					stats.UnphasedHets++
				} else {
					site.Status = SitePhased
					site.BlockID = bid
					site.H1 = sc.phase[v]
					site.H2 = 1 - site.H1
					stats.PhasedHets++
				}
			}
			if lowCoverage[v] {
				stats.LowCoverageSites++
			}
			sr.Sites[v] = site
		}

		sr.Blocks = buildBlockInfo(sc, blockIDByComp, compIDs, singletonTrioBlock, in.Variants)

		// Trio annotations may have been attached before assembly; copy them.
		if annotations, ok := ctx.siteAnnotations[sc.id]; ok {
			for v, a := range annotations {
				sr.Sites[v].PaternalH1 = a.paternalH1
				sr.Sites[v].OriginResolved = a.originResolved
				sr.Sites[v].MendelViolation = a.mendelViolation
				if a.mendelViolation {
					ctx.conflict.TrioMendelViolations++
					sr.Sites[v].Status = SiteUnphased
					sr.Sites[v].BlockID = 0
					if sc.class[v] == gtHet {
						violatingHets++
					}
				}
			}
		}

		sampleResults = append(sampleResults, sr)
	}

	res.Samples = sampleResults
	ctx.conflict.sortEntries()
	res.Conflicts = ctx.conflict
	// A site is counted as correctly phased only when it is phased and is not a
	// pedigree violation. Homozygous calls are trivially resolved and excluded
	// from this heterozygous-site metric.
	stats.CorrectlyPhased = stats.PhasedHets - violatingHets
	res.Stats = stats
}

// siteAnnotation carries trio-derived marks attached by applyTrio.
type siteAnnotation struct {
	paternalH1      bool
	originResolved  bool
	mendelViolation bool
}

func buildBlockInfo(sc *sampleCtx, blockIDByComp map[int]int, compIDs []int,
	singleton map[int]int, variants []Variant) []BlockInfo {
	// Gather het-list positions per block in canonical traversal order.
	positionsByBlock := make(map[int][]int)
	for k, cid := range compIDs {
		bid := blockIDByComp[cid]
		if bid == 0 {
			if sb, isSingleton := singleton[k]; isSingleton {
				bid = sb
			} else {
				continue
			}
		}
		positionsByBlock[bid] = append(positionsByBlock[bid], k)
	}

	blocks := make([]BlockInfo, 0, len(positionsByBlock))
	for bid := 1; bid <= len(positionsByBlock); bid++ {
		positions, ok := positionsByBlock[bid]
		if !ok {
			continue
		}
		startV := sc.hetIdx[positions[0]]
		endV := sc.hetIdx[positions[len(positions)-1]]

		// Block confidence is the product of per-edge posterior confidences
		// along consecutive het positions inside the block.
		confidence := 1.0
		supportingReads := 0
		for x := 0; x+1 < len(positions); x++ {
			if positions[x+1] != positions[x]+1 {
				continue // component contains only accepted edges anyway
			}
			edge := sc.edges[positions[x]]
			if edge == nil {
				continue
			}
			confidence *= edge.confidence
			supportingReads += edge.votesFlip + edge.votesSame
		}
		blocks = append(blocks, BlockInfo{
			SampleID:          sc.id,
			BlockID:           bid,
			StartVariantIdx:   startV,
			EndVariantIdx:     endV,
			Hets:              len(positions),
			SupportingReads:   supportingReads,
			Confidence:        roundStable(confidence),
			ExpectedErrorRate: roundStable(1 - confidence),
		})
	}
	return blocks
}

// roundStable rounds to 12 decimal places so float rounding artifacts never
// produce differing output bytes across runs or hosts.
func roundStable(x float64) float64 {
	return math.Round(x*1e12) / 1e12
}

func validatePrior(in *Input, prior *Result) error {
	if len(prior.Variants) != len(in.Variants) {
		return fmt.Errorf("rephase: prior has %d variants, input has %d",
			len(prior.Variants), len(in.Variants))
	}
	for i := range in.Variants {
		if prior.Variants[i].ID != in.Variants[i].ID ||
			prior.Variants[i].Chrom != in.Variants[i].Chrom ||
			prior.Variants[i].Pos != in.Variants[i].Pos {
			return fmt.Errorf("rephase: variant %d mismatch between prior and input", i)
		}
	}
	sampleByID := make(map[string]Sample, len(in.Samples))
	for _, s := range in.Samples {
		sampleByID[s.ID] = s
	}
	if len(prior.Samples) != len(in.Samples) {
		return fmt.Errorf("rephase: sample count mismatch")
	}
	for _, sr := range prior.Samples {
		s, ok := sampleByID[sr.ID]
		if !ok {
			return fmt.Errorf("rephase: unknown sample %q in prior", sr.ID)
		}
		if len(sr.Sites) != len(in.Variants) {
			return fmt.Errorf("rephase: sample %q site count mismatch", sr.ID)
		}
		for v, site := range sr.Sites {
			call := s.GT[v]
			switch site.Status {
			case SitePhased:
				if site.H1 != -1 && site.H2 != -1 {
					if !allelesMatch(call, site.H1, site.H2) {
						return fmt.Errorf("rephase: sample %q variant %d phased alleles incompatible with genotype",
							sr.ID, v)
					}
				}
			case SiteUnphased, SiteMissing:
				// allowed: rephase only needs identity of already-phased entries
			default:
				return fmt.Errorf("rephase: unknown site status %q", site.Status)
			}
		}
	}
	return nil
}

func reconstructFromPrior(in *Input, prior *Result) *Result {
	// Re-run the pipeline from the identical input and config, then overlay the
	// prior phasing assignment on heterozygous sites so that rephasing a service
	// result yields byte-identical haplotypes and blocks.
	ctx, _ := buildContext(in)
	phaseSamples(ctx)
	if in.Trio != nil {
		applyTrio(ctx)
	}
	assembleResult(ctx)

	priorSites := make(map[string][]SitePhase, len(prior.Samples))
	priorBlocks := make(map[string][]BlockInfo, len(prior.Samples))
	for _, sr := range prior.Samples {
		priorSites[sr.ID] = sr.Sites
		priorBlocks[sr.ID] = sr.Blocks
	}
	for i := range ctx.result.Samples {
		sr := &ctx.result.Samples[i]
		ps, ok := priorSites[sr.ID]
		if !ok {
			continue
		}
		for v := range sr.Sites {
			if ps[v].Status == SitePhased && ps[v].H1 >= 0 {
				sr.Sites[v] = ps[v]
			}
		}
		sr.Blocks = priorBlocks[sr.ID]
	}
	// A rephase re-emits the same run: preserve the prior aggregates so the
	// output is byte identical, not merely haplotype-equivalent.
	ctx.result.Stats = prior.Stats
	ctx.result.Conflicts = prior.Conflicts
	return ctx.result
}

// allelesMatch reports whether an unordered diploid call contains the two
// phased alleles, treating -1 as a wildcard for missing data.
func allelesMatch(call GenotypeCall, h1, h2 int8) bool {
	if !call.Present {
		return true
	}
	possibilities := [][2]int8{{call.A1, call.A2}}
	if call.A1 != call.A2 {
		possibilities = append(possibilities, [2]int8{call.A2, call.A1})
	}
	for _, p := range possibilities {
		if alleleEqual(p[0], h1) && alleleEqual(p[1], h2) {
			return true
		}
	}
	return false
}

func alleleEqual(callAllele, phasedAllele int8) bool {
	return callAllele == -1 || callAllele == phasedAllele
}

func (sc *sampleCtx) sampleID() string {
	return sc.id
}

func (ctx *engineCtx) addConflict(e ConflictEntry) {
	ctx.conflict.Entries = append(ctx.conflict.Entries, e)
}

func max0(n int) int {
	if n < 0 {
		return 0
	}
	return n
}

// sortEntries orders conflict entries deterministically for stable output.
func (r *ConflictReport) sortEntries() {
	sort.Slice(r.Entries, func(i, j int) bool {
		a, b := r.Entries[i], r.Entries[j]
		if a.Kind != b.Kind {
			return a.Kind < b.Kind
		}
		if a.SampleID != b.SampleID {
			return a.SampleID < b.SampleID
		}
		if a.VariantI != b.VariantI {
			return a.VariantI < b.VariantI
		}
		if a.VariantJ != b.VariantJ {
			return a.VariantJ < b.VariantJ
		}
		return a.ReadID < b.ReadID
	})
}

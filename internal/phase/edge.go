package phase

import (
	"sort"

	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

// siteInfo caches the per-sample classification of one variant site.
type siteInfo struct {
	index    int
	hetero   bool
	present  bool
	ploidy   int
	violated bool
	a0       string
	a1       string
	depth    int
	lowCov   bool
}

// voteTally holds read- and trio-derived votes between adjacent het sites.
// same = both sites carry the canonical A0 (or both A1) on H1; opp = flipped.
type voteTally struct {
	fromHet  int // ordinal of the first het site in the pair
	readSame int
	readOpp  int
	trioSame int
	trioOpp  int
}

func (v voteTally) totals() (same, opp int) {
	return v.readSame + v.trioSame, v.readOpp + v.trioOpp
}

// classifySites computes canonical alleles, hetero status and read depth.
func classifySites(ds *model.Dataset, sample string, minDepth int) []siteInfo {
	gts := ds.GT[sample]
	sites := make([]siteInfo, len(ds.Variants))
	depths := make([]int, len(ds.Variants))
	for _, read := range ds.Reads {
		if read.Sample != sample {
			continue
		}
		for _, obs := range read.Observations {
			depths[obs.VariantIdx]++
		}
	}
	for i, variant := range ds.Variants {
		gt := gts[i]
		info := siteInfo{
			index:   i,
			present: gt.Present,
			ploidy:  gt.Ploidy,
			a0:      variant.Ref,
			a1:      variant.Alt,
			depth:   depths[i],
			lowCov:  depths[i] < minDepth,
		}
		if gt.Present && gt.Ploidy == 2 && gt.A1 != gt.A2 {
			info.hetero = true
			info.a0, info.a1 = canonicalAlleles(variant, gt.A1, gt.A2)
		}
		sites[i] = info
	}
	return sites
}

// canonicalAlleles returns the two distinct alleles with Ref first as A0.
func canonicalAlleles(variant model.Variant, x, y string) (string, string) {
	a0, a1 := x, y
	if a0 > a1 {
		a0, a1 = a1, a0
	}
	if variant.Ref == a1 {
		a0, a1 = a1, a0
	}
	return a0, a1
}

// collectReadEdges aggregates same/opposite read votes between adjacent
// heterozygous sites. Only a pair of consecutive het sites within a read
// forms one direct edge vote.
func collectReadEdges(ds *model.Dataset, sample string, sites []siteInfo) []voteTally {
	hetIdx := map[int]int{}
	ordinal := 0
	for _, info := range sites {
		if info.hetero && !info.violated {
			hetIdx[info.index] = ordinal
			ordinal++
		}
	}
	edges := make([]voteTally, 0)
	edgeAt := map[int]int{}

	alleleBit := func(info siteInfo, allele string) int {
		if allele == info.a1 {
			return 1
		}
		return 0
	}

	for _, read := range ds.Reads {
		if read.Sample != sample {
			continue
		}
		type obs struct {
			siteOrdinal int
			bit         int
		}
		var seq []obs
		for _, o := range read.Observations {
			h, ok := hetIdx[o.VariantIdx]
			if !ok {
				continue
			}
			seq = append(seq, obs{siteOrdinal: h, bit: alleleBit(sites[o.VariantIdx], o.Allele)})
		}
		sort.SliceStable(seq, func(i, j int) bool { return seq[i].siteOrdinal < seq[j].siteOrdinal })
		for k := 1; k < len(seq); k++ {
			if seq[k].siteOrdinal != seq[k-1].siteOrdinal+1 {
				continue
			}
			h := seq[k-1].siteOrdinal
			pos, ok := edgeAt[h]
			if !ok {
				pos = len(edges)
				edgeAt[h] = pos
				edges = append(edges, voteTally{fromHet: h})
			}
			if seq[k-1].bit == seq[k].bit {
				edges[pos].readSame++
			} else {
				edges[pos].readOpp++
			}
		}
	}
	return edges
}

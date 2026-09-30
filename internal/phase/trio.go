package phase

import "github.com/guoyyit/haplotype-phasing-go/internal/model"

// trioEdgeWeight is the effective vote weight of one Mendelian-inheritance
// edge relative to one read-pair vote. Trio evidence is high confidence.
const trioEdgeWeight = 3

// alleleSet returns the set of observed alleles at one site, plus a flag
// for a completely missing genotype.
func alleleSet(gt model.Genotype) (map[string]bool, bool) {
	if !gt.Present || gt.Ploidy == 0 {
		return map[string]bool{}, true
	}
	set := map[string]bool{gt.A1: true}
	if gt.Ploidy == 2 {
		set[gt.A2] = true
	}
	return set, false
}

// mendelianViolation reports whether a diploid child genotype is impossible
// given the two parental genotypes. Each of the child's two alleles must be
// traceable to a distinct parent (one paternal, one maternal), and the
// parental transmission must respect homozygote dosage. Missing data on any
// side is uninformative (returns false).
func mendelianViolation(child, father, mother model.Genotype) bool {
	if !child.Present || child.Ploidy < 2 {
		return false
	}
	fSet, fMiss := alleleSet(father)
	mSet, mMiss := alleleSet(mother)
	if fMiss || mMiss {
		return false
	}
	// Enumerate the two ordered child allele assignments (A1,A2) and (A2,A1):
	// first allele paternal, second maternal.
	assignments := [2][2]string{{child.A1, child.A2}, {child.A2, child.A1}}
	for _, as := range assignments {
		if fSet[as[0]] && mSet[as[1]] {
			return false
		}
	}
	return true
}

// paternalBit returns, for one heterozygous child site, the canonical bit
// (0=A0 paternal, 1=A1 paternal) when Mendelian inheritance forces it.
// ok=false means the inheritance is uninformative (e.g. both parents het).
func paternalBit(child, father, mother model.Genotype, a0, a1 string) (bit int, ok bool) {
	if !child.Present || child.Ploidy < 2 {
		return 0, false
	}
	fSet, fMiss := alleleSet(father)
	mSet, mMiss := alleleSet(mother)
	if fMiss || mMiss {
		return 0, false
	}
	possible := map[string]bool{}
	for _, allele := range []string{a0, a1} {
		other := a0
		if allele == a0 {
			other = a1
		}
		if fSet[allele] && mSet[other] {
			possible[allele] = true
		}
	}
	// Only when exactly one allele can be paternal is the phase forced.
	if possible[a0] && !possible[a1] {
		return 0, true
	}
	if possible[a1] && !possible[a0] {
		return 1, true
	}
	return 0, false
}

// trioEvidence marks per-site violations and returns, indexed by het-site
// ordinal, the forced paternal bit (valid only where forced is true).
func trioEvidence(ds *model.Dataset, trio model.Trio, sites []siteInfo) (forced []bool, bits []int) {
	n := 0
	for _, info := range sites {
		if info.hetero && !info.violated {
			n++
		}
	}
	forced = make([]bool, n)
	bits = make([]int, n)

	fatherGT := ds.GT[trio.Father]
	motherGT := ds.GT[trio.Mother]
	childGT := ds.GT[trio.Child]
	ordinal := 0
	for i := range sites {
		fgt, mgt, cgt := gtAt(fatherGT, i), gtAt(motherGT, i), gtAt(childGT, i)
		// Violation check applies to every diploid child site.
		if mendelianViolation(cgt, fgt, mgt) {
			sites[i].violated = true
		}
		info := sites[i]
		if !info.hetero || info.violated {
			continue
		}
		if bit, ok := paternalBit(cgt, fgt, mgt, info.a0, info.a1); ok {
			forced[ordinal] = true
			bits[ordinal] = bit
		}
		ordinal++
	}
	return forced, bits
}

// gtAt returns the genotype at i, or a fully-missing genotype when the
// sample has no (or a shorter) genotype vector.
func gtAt(gts []model.Genotype, i int) model.Genotype {
	if i < 0 || i >= len(gts) {
		return model.Genotype{}
	}
	return gts[i]
}

// addTrioEdges casts one weighted trio vote on every consecutive pair of
// heterozygous child sites whose paternal bits are both forced.
func addTrioEdges(edges []voteTally, forced []bool, bits []int) []voteTally {
	edgeAt := map[int]int{}
	for i := range edges {
		edgeAt[edges[i].fromHet] = i
	}
	for h := 0; h+1 < len(forced); h++ {
		if !forced[h] || !forced[h+1] {
			continue
		}
		pos, ok := edgeAt[h]
		if !ok {
			pos = len(edges)
			edgeAt[h] = pos
			edges = append(edges, voteTally{fromHet: h})
		}
		if bits[h] == bits[h+1] {
			edges[pos].trioSame += trioEdgeWeight
		} else {
			edges[pos].trioOpp += trioEdgeWeight
		}
	}
	return edges
}

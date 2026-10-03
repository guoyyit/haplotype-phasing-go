package phase

import (
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

func trioDataset() *model.Dataset {
	return &model.Dataset{
		Variants: []model.Variant{
			{ID: "v1", Ref: "A", Alt: "G"},
			{ID: "v2", Ref: "C", Alt: "T"},
		},
		Samples: []string{"child", "father", "mother"},
		GT: map[string][]model.Genotype{
			"child": {
				{Present: true, Ploidy: 2, A1: "A", A2: "G"},
				{Present: true, Ploidy: 2, A1: "C", A2: "T"},
			},
			"father": {
				{Present: true, Ploidy: 2, A1: "A", A2: "A"},
				{Present: true, Ploidy: 2, A1: "C", A2: "T"},
			},
			"mother": {
				{Present: true, Ploidy: 2, A1: "G", A2: "G"},
				{Present: true, Ploidy: 2, A1: "C", A2: "T"},
			},
		},
		Trios: []model.Trio{{Child: "child", Father: "father", Mother: "mother"}},
	}
}

// Child H1 must carry the paternal allele when inheritance forces it.
func TestChildHaplotypeInheritsPaternal(t *testing.T) {
	ds := trioDataset()
	res := Run(ds, config.Default())
	child := res.Samples["child"]
	// father A/A, mother G/G => paternal allele at v1 is A (= A0).
	if child.Sites[0].H1 != model.AlleleA0 || child.Sites[0].H2 != model.AlleleA1 {
		t.Fatalf("expected H1=paternal A0 / H2=maternal A1, got %v/%v", child.Sites[0].H1, child.Sites[0].H2)
	}
	if child.ViolationHets != 0 {
		t.Fatalf("expected no violations, got %d", child.ViolationHets)
	}
	if res.Stats.TrioViolations != 0 {
		t.Fatalf("expected zero trio violations, got %d", res.Stats.TrioViolations)
	}
}

// A genotype impossible under Mendelian inheritance is counted and marked,
// and must not be counted as correctly phased.
func TestMendelianViolationMarkedAndUnphased(t *testing.T) {
	ds := &model.Dataset{
		Variants: []model.Variant{{ID: "v1", Ref: "A", Alt: "G"}},
		Samples:  []string{"c", "f", "m"},
		GT: map[string][]model.Genotype{
			"c": {{Present: true, Ploidy: 2, A1: "A", A2: "A"}},
			"f": {{Present: true, Ploidy: 2, A1: "G", A2: "G"}},
			"m": {{Present: true, Ploidy: 2, A1: "G", A2: "G"}},
		},
		Trios: []model.Trio{{Child: "c", Father: "f", Mother: "m"}},
	}
	res := Run(ds, config.Default())
	site := res.Samples["c"].Sites[0]
	if !site.Violation {
		t.Fatal("expected the impossible child genotype to be flagged")
	}
	if site.Phased {
		t.Fatal("violating site must not be phased")
	}
	if res.Stats.TrioViolations != 1 {
		t.Fatalf("expected 1 trio violation, got %d", res.Stats.TrioViolations)
	}
	if res.Samples["c"].PhasedHets != 0 {
		t.Fatal("violating site must not be counted as phased")
	}
}

// Missing parental data never produces a violation (uninformative).
func TestMissingParentsAreUninformative(t *testing.T) {
	ds := trioDataset()
	ds.GT["father"][0] = model.Genotype{Present: true, Ploidy: 0}
	res := Run(ds, config.Default())
	if res.Stats.TrioViolations != 0 {
		t.Fatalf("missing parent must not trigger violation, got %d", res.Stats.TrioViolations)
	}
}

// A violated heterozygous site must record the canonical allele calls so the
// genotype can be reconstructed on re-read; the phase stays unphased.
func TestViolatedHetRecordsCanonicalAlleleCalls(t *testing.T) {
	ds := &model.Dataset{
		Variants: []model.Variant{{ID: "v1", Ref: "A", Alt: "G"}},
		Samples:  []string{"c", "f", "m"},
		GT: map[string][]model.Genotype{
			"c": {{Present: true, Ploidy: 2, A1: "A", A2: "G"}},
			"f": {{Present: true, Ploidy: 2, A1: "G", A2: "G"}},
			"m": {{Present: true, Ploidy: 2, A1: "G", A2: "G"}},
		},
		Trios: []model.Trio{{Child: "c", Father: "f", Mother: "m"}},
	}
	res := Run(ds, config.Default())
	site := res.Samples["c"].Sites[0]
	if !site.Violation || !site.Hetero || site.Phased {
		t.Fatalf("violated het must be flagged, hetero and unphased: %+v", site)
	}
	if site.H1 != model.AlleleA0 || site.H2 != model.AlleleA1 {
		t.Fatalf("violated het must record canonical allele calls A0/A1: %v/%v", site.H1, site.H2)
	}
	if res.Stats.TrioViolations != 1 {
		t.Fatalf("expected 1 trio violation, got %d", res.Stats.TrioViolations)
	}
}

// A violated homozygous site must record the homozygous allele calls.
func TestViolatedHomoRecordsAlleleCalls(t *testing.T) {
	ds := &model.Dataset{
		Variants: []model.Variant{{ID: "v1", Ref: "C", Alt: "T"}},
		Samples:  []string{"c", "f", "m"},
		GT: map[string][]model.Genotype{
			"c": {{Present: true, Ploidy: 2, A1: "C", A2: "C"}},
			"f": {{Present: true, Ploidy: 2, A1: "T", A2: "T"}},
			"m": {{Present: true, Ploidy: 2, A1: "T", A2: "T"}},
		},
		Trios: []model.Trio{{Child: "c", Father: "f", Mother: "m"}},
	}
	res := Run(ds, config.Default())
	site := res.Samples["c"].Sites[0]
	if !site.Violation || site.Hetero || site.Phased {
		t.Fatalf("violated homo must be flagged, not hetero and unphased: %+v", site)
	}
	if site.H1 != model.AlleleA0 || site.H2 != model.AlleleA0 {
		t.Fatalf("violated homo must record homozygous allele calls: %v/%v", site.H1, site.H2)
	}
}

// Weighted trio evidence overrides a single contradicting read pair and the
// read-vs-trio orientation disagreement is reported, not silently dropped.
func TestTrioEdgeOverridesSingleReadAndCountsConflict(t *testing.T) {
	ds := &model.Dataset{
		Variants: []model.Variant{
			{ID: "v1", Ref: "A", Alt: "G"},
			{ID: "v2", Ref: "C", Alt: "T"},
		},
		Samples: []string{"child", "father", "mother"},
		GT: map[string][]model.Genotype{
			// Both sites: father homozygous for one allele, mother for the
			// other, so the paternal bit is uniquely forced at both sites.
			"child": {
				{Present: true, Ploidy: 2, A1: "A", A2: "G"},
				{Present: true, Ploidy: 2, A1: "C", A2: "T"},
			},
			"father": {
				{Present: true, Ploidy: 2, A1: "A", A2: "A"},
				{Present: true, Ploidy: 2, A1: "T", A2: "T"},
			},
			"mother": {
				{Present: true, Ploidy: 2, A1: "G", A2: "G"},
				{Present: true, Ploidy: 2, A1: "C", A2: "C"},
			},
		},
		Trios: []model.Trio{{Child: "child", Father: "father", Mother: "mother"}},
	}
	// v1 paternal=A(bit0), v2 paternal=T(bit1): trio says "opp".
	// Read says "same" (A then C=bit0).
	ds.Reads = []model.Read{
		{ID: "r1", Sample: "child", Observations: []model.ReadObservation{
			{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "C"},
		}},
	}
	cfg := config.Default()
	cfg.MinEdgeSupport = 1 // only one read, but trio weight supplies strength
	res := Run(ds, cfg)
	if res.Stats.ReadEdgeConflicts == 0 {
		t.Fatal("read-vs-trio disagreement must be counted")
	}
	child := res.Samples["child"]
	// Trio wins (weight 3 vs 1 read): H1 paternal at v1 = A0, v2 = A1.
	if child.Sites[0].H1 != model.AlleleA0 || child.Sites[1].H1 != model.AlleleA1 {
		t.Fatalf("trio majority must orient H1 paternal: %v %v",
			child.Sites[0].H1, child.Sites[1].H1)
	}
}

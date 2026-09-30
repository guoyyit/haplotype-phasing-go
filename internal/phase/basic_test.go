package phase

import (
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

// Without reads or trios, every het is a deterministic SINGLETON phased
// with A0 placed on H1; output is still complete.
func TestGenotypeOnlySingletonPhasing(t *testing.T) {
	ds := &model.Dataset{
		Variants: []model.Variant{
			{ID: "v1", Ref: "A", Alt: "G"},
			{ID: "v2", Ref: "C", Alt: "T"},
		},
		Samples: []string{"s"},
		GT: map[string][]model.Genotype{
			"s": {
				{Present: true, Ploidy: 2, A1: "A", A2: "G"},
				{Present: true, Ploidy: 2, A1: "T", A2: "C"},
			},
		},
	}
	res := Run(ds, config.Default())
	sr := res.Samples["s"]
	if len(sr.Blocks) != 2 {
		t.Fatalf("expected 2 singleton blocks, got %d", len(sr.Blocks))
	}
	if sr.PhasedHets != 2 || sr.UnphasedHets != 0 {
		t.Fatalf("both hets should be phased, got %d/%d", sr.PhasedHets, sr.UnphasedHets)
	}
	// Canonical default orientation: A0 on H1 for every singleton.
	if sr.Sites[0].H1 != model.AlleleA0 || sr.Sites[1].H1 != model.AlleleA0 {
		t.Fatalf("default singleton orientation must be A0 on H1: %v %v",
			sr.Sites[0].H1, sr.Sites[1].H1)
	}
}

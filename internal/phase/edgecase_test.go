package phase

import (
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

// A trio referencing parents without genotype rows must not panic; the
// missing parents are uninformative.
func TestTrioWithAbsentParentNoPanic(t *testing.T) {
	ds := &model.Dataset{
		Variants: []model.Variant{{ID: "v1", Ref: "A", Alt: "G"}},
		Samples:  []string{"child"},
		GT: map[string][]model.Genotype{
			"child": {{Present: true, Ploidy: 2, A1: "A", A2: "G"}},
		},
		Trios: []model.Trio{{Child: "child", Father: "ghost-dad", Mother: "ghost-mom"}},
	}
	res := Run(ds, config.Default())
	if res.Samples["child"] == nil {
		t.Fatal("child result missing")
	}
}

// A read referencing a sample absent from genotypes must not panic.
func TestReadForUnknownSampleNoPanic(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("phasing panicked: %v", r)
		}
	}()
	ds := &model.Dataset{
		Variants: []model.Variant{{ID: "v1", Ref: "A", Alt: "G"}, {ID: "v2", Ref: "C", Alt: "T"}},
		Samples:  []string{"s"},
		GT: map[string][]model.Genotype{
			"s": {
				{Present: true, Ploidy: 2, A1: "A", A2: "G"},
				{Present: true, Ploidy: 2, A1: "C", A2: "T"},
			},
		},
		Reads: []model.Read{
			{ID: "r", Sample: "ghost", Observations: []model.ReadObservation{
				{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "C"},
			}},
		},
	}
	Run(ds, config.Default())
}

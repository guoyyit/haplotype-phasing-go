package phase

import (
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

// Same input and config must produce identical results repeatedly.
func TestDeterministicRepeatedRuns(t *testing.T) {
	ds := heteroDataset()
	ds.Reads = []model.Read{
		{ID: "r1", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "T"}}},
		{ID: "r2", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "C"}}},
	}
	cfg := config.Default()
	first := Run(ds, cfg)
	for i := 0; i < 5; i++ {
		again := Run(ds, cfg)
		if !sameResultShape(first, again) {
			t.Fatalf("run %d differs from the first run", i+1)
		}
	}
}

// A fixed seed must reproduce the exact tie-break; a different seed is
// allowed to flip the tie but must itself be reproducible.
func TestSeedReproducibility(t *testing.T) {
	ds := heteroDataset()
	ds.Reads = []model.Read{
		{ID: "r1", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "T"}}},
		{ID: "r2", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "G"}, {VariantIdx: 1, Allele: "T"}}},
	}

	seed1 := config.Default()
	seed1.Seed = 42
	first := Run(ds, seed1).Samples["s"].Sites[1].H1
	for i := 0; i < 3; i++ {
		if got := Run(ds, seed1).Samples["s"].Sites[1].H1; got != first {
			t.Fatalf("seed 42 must reproduce, run %d got %v want %v", i, got, first)
		}
	}

	seed2 := config.Default()
	seed2.Seed = 99
	a := Run(ds, seed2).Samples["s"].Sites[1].H1
	b := Run(ds, seed2).Samples["s"].Sites[1].H1
	if a != b {
		t.Fatal("same non-default seed must be reproducible")
	}
}

func sameResultShape(a, b *model.Result) bool {
	if len(a.Order) != len(b.Order) || a.Stats != b.Stats {
		return false
	}
	for _, sample := range a.Order {
		sa, sb := a.Samples[sample], b.Samples[sample]
		if len(sa.Sites) != len(sb.Sites) || len(sa.Blocks) != len(sb.Blocks) {
			return false
		}
		for i := range sa.Sites {
			x, y := sa.Sites[i], sb.Sites[i]
			if x.H1 != y.H1 || x.H2 != y.H2 || x.Phased != y.Phased || x.Violation != y.Violation {
				return false
			}
		}
	}
	return true
}

// Adding another sample must not change an existing sample's tie-breaks.
func TestSampleSetIndependence(t *testing.T) {
	mk := func(samples ...string) *model.Dataset {
		ds := &model.Dataset{
			Variants: []model.Variant{{ID: "v1", Ref: "A", Alt: "G"}, {ID: "v2", Ref: "C", Alt: "T"}},
			GT:       map[string][]model.Genotype{},
		}
		for _, name := range samples {
			ds.Samples = append(ds.Samples, name)
			ds.GT[name] = []model.Genotype{
				{Present: true, Ploidy: 2, A1: "A", A2: "G"},
				{Present: true, Ploidy: 2, A1: "C", A2: "T"},
			}
			ds.Reads = append(ds.Reads,
				model.Read{ID: name + "-r1", Sample: name, Observations: []model.ReadObservation{
					{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "T"},
				}},
				model.Read{ID: name + "-r2", Sample: name, Observations: []model.ReadObservation{
					{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "C"},
				}},
			)
		}
		return ds
	}
	cfg := config.Default()
	one := Run(mk("alpha"), cfg).Samples["alpha"].Sites[1].H1
	two := Run(mk("alpha", "beta"), cfg).Samples["alpha"].Sites[1].H1
	if one != two {
		t.Fatalf("alpha tie-break changed when beta was added: %v vs %v", one, two)
	}
}

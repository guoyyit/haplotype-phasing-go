package phase

import (
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

func heteroDataset() *model.Dataset {
	return &model.Dataset{
		Variants: []model.Variant{
			{ID: "v1", Ref: "A", Alt: "G"},
			{ID: "v2", Ref: "C", Alt: "T"},
		},
		Samples: []string{"s"},
		GT: map[string][]model.Genotype{
			"s": {
				{Present: true, Ploidy: 2, A1: "A", A2: "G"},
				{Present: true, Ploidy: 2, A1: "C", A2: "T"},
			},
		},
	}
}

// Unanimous "opposite" reads connect a block and flip H1 at the second site.
func TestReadEvidenceOppositeFlipsHaplotype(t *testing.T) {
	ds := heteroDataset()
	ds.Reads = []model.Read{
		{ID: "r1", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "T"}}},
		{ID: "r2", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "G"}, {VariantIdx: 1, Allele: "C"}}},
	}
	res := Run(ds, config.Default())
	sr := res.Samples["s"]
	if len(sr.Blocks) != 1 {
		t.Fatalf("expected 1 connected block, got %d", len(sr.Blocks))
	}
	if sr.Sites[1].H1 != model.AlleleA1 {
		t.Fatalf("expected flipped H1=A1 at site 2, got %v", sr.Sites[1].H1)
	}
	if res.Stats.ReadEdgeConflicts != 0 {
		t.Fatalf("expected no conflicts, got %d", res.Stats.ReadEdgeConflicts)
	}
}

// Conflicting read votes (majority vs minority) are counted, not dropped.
func TestReadConflictCountedAndAdjudicated(t *testing.T) {
	ds := heteroDataset()
	ds.Reads = []model.Read{
		{ID: "r1", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "T"}}},
		{ID: "r2", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "G"}, {VariantIdx: 1, Allele: "C"}}},
		{ID: "r3", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "C"}}},
	}
	res := Run(ds, config.Default())
	if res.Stats.ReadConflictsAdjudicated == 0 {
		t.Fatal("expected at least one adjudicated conflict")
	}
	if res.Stats.ReadEdgeConflicts == 0 {
		t.Fatal("expected conflicts to be counted")
	}
	// 2 opp vs 1 same: strength 2/3 < 0.8 default -> edge must break.
	if res.Stats.BreaksByThreshold != 1 {
		t.Fatalf("expected 1 threshold break, got %d", res.Stats.BreaksByThreshold)
	}
	if got := len(res.Samples["s"].Blocks); got != 2 {
		t.Fatalf("expected 2 blocks after conflict break, got %d", got)
	}
}

// Exact ties use the seeded tie-break deterministically and are reported.
func TestTieIsDeterministicAndReported(t *testing.T) {
	ds := heteroDataset()
	ds.Reads = []model.Read{
		{ID: "r1", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "T"}}},
		{ID: "r2", Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "C"}}},
	}
	cfg := config.Default()
	first := Run(ds, cfg)
	second := Run(ds, cfg)
	if first.Samples["s"].Sites[1].H1 != second.Samples["s"].Sites[1].H1 {
		t.Fatal("same seed must reproduce the tie-break")
	}
	if first.Stats.ReadConflictsAdjudicated == 0 {
		t.Fatal("the tie must be reported as adjudicated")
	}
}

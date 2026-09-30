package phase

import (
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

func threeHetDataset(votesPerEdge int) *model.Dataset {
	ds := &model.Dataset{
		Variants: []model.Variant{
			{ID: "v1", Ref: "A", Alt: "G"},
			{ID: "v2", Ref: "C", Alt: "T"},
			{ID: "v3", Ref: "G", Alt: "A"},
		},
		Samples: []string{"s"},
		GT: map[string][]model.Genotype{
			"s": {
				{Present: true, Ploidy: 2, A1: "A", A2: "G"},
				{Present: true, Ploidy: 2, A1: "C", A2: "T"},
				{Present: true, Ploidy: 2, A1: "G", A2: "A"},
			},
		},
	}
	for r := 0; r < votesPerEdge; r++ {
		ds.Reads = append(ds.Reads,
			model.Read{ID: readID("a", r), Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 0, Allele: "A"}, {VariantIdx: 1, Allele: "C"}}},
			model.Read{ID: readID("b", r), Sample: "s", Observations: []model.ReadObservation{{VariantIdx: 1, Allele: "C"}, {VariantIdx: 2, Allele: "G"}}},
		)
	}
	return ds
}

func readID(prefix string, r int) string {
	return prefix + "-" + string(rune('0'+r))
}

// Strong, threshold-passing evidence keeps all three hets in one block.
func TestStrongEvidenceOneBlock(t *testing.T) {
	ds := threeHetDataset(3)
	res := Run(ds, config.Default())
	if got := len(res.Samples["s"].Blocks); got != 1 {
		t.Fatalf("expected 1 block, got %d", got)
	}
	b := res.Samples["s"].Blocks[0]
	if b.ExpectedErrorRate < 0 || b.ExpectedErrorRate > 1 {
		t.Fatalf("block error rate out of range: %v", b.ExpectedErrorRate)
	}
	if b.Confidence+b.ExpectedErrorRate < 0.999 || b.Confidence+b.ExpectedErrorRate > 1.001 {
		t.Fatalf("confidence + error must be ~1, got %v+%v", b.Confidence, b.ExpectedErrorRate)
	}
}

// Raising the strength threshold deterministically splits blocks.
func TestThresholdSplitsBlocks(t *testing.T) {
	ds := threeHetDataset(1) // single unanimous vote per edge

	loose := config.Default()
	loose.MinEdgeSupport = 1
	loose.MinEdgeStrength = 0.8
	if got := len(Run(ds, loose).Samples["s"].Blocks); got != 1 {
		t.Fatalf("loose threshold: expected 1 block, got %d", got)
	}

	strict := config.Default()
	strict.MinEdgeSupport = 5
	res := Run(ds, strict)
	if got := len(res.Samples["s"].Blocks); got != 3 {
		t.Fatalf("strict support threshold: expected 3 singleton blocks, got %d", got)
	}
	if res.Stats.BreaksByThreshold != 2 {
		t.Fatalf("expected 2 threshold breaks, got %d", res.Stats.BreaksByThreshold)
	}
}

// Isolated het with no evidence is a reported SINGLETON, not dropped.
func TestSingletonBlockReported(t *testing.T) {
	ds := heteroDataset() // no reads
	res := Run(ds, config.Default())
	blocks := res.Samples["s"].Blocks
	if len(blocks) != 2 {
		t.Fatalf("expected 2 singleton blocks, got %d", len(blocks))
	}
	for _, b := range blocks {
		if b.Status != model.BlockStatusSingleton {
			t.Fatalf("expected SINGLETON status, got %s", b.Status)
		}
	}
}

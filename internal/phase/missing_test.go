package phase

import (
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

// Completely and partially missing genotypes must not abort processing and
// must be reported as MISSING / PARTIAL with empty or single allele calls.
func TestMissingAndPartialGenotypesDoNotFail(t *testing.T) {
	ds := &model.Dataset{
		Variants: []model.Variant{
			{ID: "v1", Ref: "A", Alt: "G"},
			{ID: "v2", Ref: "C", Alt: "T"},
			{ID: "v3", Ref: "G", Alt: "A"},
		},
		Samples: []string{"s"},
		GT: map[string][]model.Genotype{
			"s": {
				{Present: true, Ploidy: 0},
				{Present: true, Ploidy: 1, A1: "C"},
				{Present: true, Ploidy: 2, A1: "G", A2: "A"},
			},
		},
	}
	res := Run(ds, config.Default())
	sr := res.Samples["s"]
	if sr.Sites[0].H1 != model.AlleleMissing || sr.Sites[0].H2 != model.AlleleMissing {
		t.Fatal("fully missing site must have no allele calls")
	}
	if sr.Sites[1].H1 == model.AlleleMissing || sr.Sites[1].H2 != model.AlleleMissing {
		t.Fatal("partial site must carry exactly one allele call")
	}
	if !sr.Sites[2].Phased {
		t.Fatal("the valid het site should still be phased")
	}
	if res.Stats.MissingGenotypes < 1 {
		t.Fatal("missing genotype must be counted")
	}
}

// Low-coverage sites are flagged but still phasable.
func TestLowCoverageFlaggedButPhasable(t *testing.T) {
	ds := heteroDataset()
	cfg := config.Default()
	cfg.MinSiteDepth = 5 // no read reaches depth 5
	res := Run(ds, cfg)
	if res.Stats.LowCoverageSites != 2 {
		t.Fatalf("expected 2 low-coverage sites, got %d", res.Stats.LowCoverageSites)
	}
	if !res.Samples["s"].Sites[0].LowCov {
		t.Fatal("site should be flagged low coverage")
	}
}

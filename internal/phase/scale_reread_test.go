package phase_test

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
	"github.com/guoyyit/haplotype-phasing-go/internal/output"
	"github.com/guoyyit/haplotype-phasing-go/internal/parseinput"
	"github.com/guoyyit/haplotype-phasing-go/internal/phase"
)

// TestScaleBoundaryReread exercises the full write/re-read/re-phase/write
// verification path at two size tiers and asserts:
//   - the boundary run finishes well inside the documented time budget;
//   - the regenerated file is byte-identical;
//   - verification wall time grows smoothly: quadrupling the work must not
//     take more than 10x the time. The slack absorbs GC/scheduler jitter on
//     small machines while still ruling out any quadratic blow-up.
func TestScaleBoundaryReread(t *testing.T) {
	if testing.Short() {
		t.Skip("scale boundary reread skipped in -short mode")
	}
	tier := func(nVariants, nSamples, nReads int) (time.Duration, int) {
		ds := makeExternalScaleDataset(nVariants, nSamples, nReads)
		cfg := config.Default()
		res := phase.Run(ds, cfg)

		start := time.Now()
		var first bytes.Buffer
		if err := output.WritePhased(&first, ds, res, cfg); err != nil {
			t.Fatal(err)
		}
		path := filepath.Join(t.TempDir(), "scaled.tsv")
		if err := os.WriteFile(path, first.Bytes(), 0644); err != nil {
			t.Fatal(err)
		}
		ds2, res2, cfg2, _, version, err := parseinput.LoadPhased(path)
		if err != nil {
			t.Fatal(err)
		}
		var second bytes.Buffer
		if err := output.WritePhasedVersion(&second, ds2, phase.Run(ds2, cfg2), cfg2, version); err != nil {
			t.Fatal(err)
		}
		if !bytes.Equal(first.Bytes(), second.Bytes()) {
			t.Fatalf("reread mismatch at scale %d/%d", nVariants, nSamples)
		}
		if res2.Stats != res.Stats {
			t.Fatalf("stats changed across scale reread: %+v vs %+v", res.Stats, res2.Stats)
		}
		return time.Since(start), first.Len()
	}

	d1, b1 := tier(2500, 20, 25000)
	largeDur, largeBytes := tier(10000, 50, 100000)

	if largeDur > 10*time.Second {
		t.Fatalf("boundary reread took %v, exceeds documented 10s phasing budget", largeDur)
	}
	if b1 >= largeBytes {
		t.Fatal("output size should grow with scale")
	}
	// The boundary tier is 4x the work of the small tier; wall time stays
	// below 10x (near-linear growth with jitter slack, never quadratic).
	if largeDur > 10*d1 {
		t.Fatalf("verification time does not scale smoothly: %v -> %v for a 4x workload", d1, largeDur)
	}
	t.Logf("reread durations small=%v boundary=%v; bytes small=%d boundary=%d",
		d1, largeDur, b1, largeBytes)
}

// makeExternalScaleDataset mirrors the internal synthetic generator without
// reaching into package-private helpers.
func makeExternalScaleDataset(nVariants, nSamples, nReads int) *model.Dataset {
	ds := &model.Dataset{GT: map[string][]model.Genotype{}}
	for i := 0; i < nVariants; i++ {
		ds.Variants = append(ds.Variants, model.Variant{
			ID: fmt.Sprintf("v%06d", i), Chrom: "chr1", Pos: i + 1, Ref: "A", Alt: "G",
		})
	}
	for s := 0; s < nSamples; s++ {
		name := fmt.Sprintf("sample%03d", s)
		ds.Samples = append(ds.Samples, name)
		gts := make([]model.Genotype, nVariants)
		for i := range gts {
			gts[i] = model.Genotype{Present: true, Ploidy: 2, A1: "A", A2: "G"}
		}
		ds.GT[name] = gts
	}
	for r := 0; r < nReads; r++ {
		sample := ds.Samples[r%nSamples]
		v := (r / nSamples) % (nVariants - 1)
		first, second := "A", "A"
		if r%2 == 1 {
			first, second = "G", "G"
		}
		ds.Reads = append(ds.Reads, model.Read{
			ID:     fmt.Sprintf("r%07d", r),
			Sample: sample,
			Observations: []model.ReadObservation{
				{VariantIdx: v, Allele: first},
				{VariantIdx: v + 1, Allele: second},
			},
		})
	}
	return ds
}

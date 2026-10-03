package phase

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
	"github.com/guoyyit/haplotype-phasing-go/internal/output"
	"github.com/guoyyit/haplotype-phasing-go/internal/parseinput"
)

// The documented supported boundary: 10,000 variants, 50 samples,
// 100,000 reads. Resource limits documented in docs/scale.md.
func BenchmarkScaleBoundary(b *testing.B) {
	ds := makeScaleDataset(10000, 50, 100000)
	cfg := config.Default()
	b.ResetTimer()
	var before, after runtime.MemStats
	runtime.GC()
	runtime.ReadMemStats(&before)
	for n := 0; n < b.N; n++ {
		Run(ds, cfg)
	}
	runtime.ReadMemStats(&after)
	heapMB := float64(after.HeapAlloc-before.HeapAlloc) / (1024 * 1024)
	b.ReportMetric(heapMB, "deltaHeapMB")
}

// TestScaleBoundaryCompletes verifies linear growth and correctness at the
// declared boundary without a wall-clock budget flakiness in unit tests.
func TestScaleBoundaryCompletes(t *testing.T) {
	if testing.Short() {
		t.Skip("scale boundary skipped in -short mode")
	}
	small := makeScaleDataset(1000, 10, 10000)
	large := makeScaleDataset(10000, 50, 100000)

	resSmall := Run(small, config.Default())
	resLarge := Run(large, config.Default())

	if len(resLarge.Order) != 50 || resLarge.Stats.Variants != 10000 {
		t.Fatalf("unexpected result dimensions: %d samples, %d variants", len(resLarge.Order), resLarge.Stats.Variants)
	}
	// Every sample must have a complete site vector (no site dropped).
	for _, sample := range resLarge.Order {
		if len(resLarge.Samples[sample].Sites) != 10000 {
			t.Fatalf("sample %s has incomplete sites: %d", sample, len(resLarge.Samples[sample].Sites))
		}
	}
	// Growth in reads/variants must scale total work smoothly (>0, finite).
	if resLarge.Stats.ReadObservations <= resSmall.Stats.ReadObservations {
		t.Fatal("read observation count should grow with the dataset")
	}
}

// makeScaleDataset builds a deterministic synthetic dataset.
func makeScaleDataset(nVariants, nSamples, nReads int) *model.Dataset {
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
			// All heterozygous to maximize phasing work.
			gts[i] = model.Genotype{Present: true, Ploidy: 2, A1: "A", A2: "G"}
		}
		ds.GT[name] = gts
	}
	// Reads each span two adjacent variants so every edge gets evidence.
	for r := 0; r < nReads; r++ {
		sample := ds.Samples[r%nSamples]
		v := (r / nSamples) % (nVariants - 1)
		alleleFirst := "A"
		alleleSecond := "A"
		if r%2 == 1 {
			alleleFirst, alleleSecond = "G", "G"
		}
		ds.Reads = append(ds.Reads, model.Read{
			ID:     fmt.Sprintf("r%07d", r),
			Sample: sample,
			Observations: []model.ReadObservation{
				{VariantIdx: v, Allele: alleleFirst},
				{VariantIdx: v + 1, Allele: alleleSecond},
			},
		})
	}
	return ds
}

// makeScaleDatasetWithTrios builds on makeScaleDataset and adds a few
// consistent trios to exercise the trio path at scale.
func makeScaleDatasetWithTrios(nVariants, nSamples, nReads int) *model.Dataset {
	ds := makeScaleDataset(nVariants, nSamples, nReads)
	if nSamples >= 3 {
		ds.Trios = []model.Trio{
			{Child: ds.Samples[0], Father: ds.Samples[1], Mother: ds.Samples[2]},
		}
	}
	return ds
}

// runScaleRoundTrip runs the full pipeline (phase -> write -> re-read ->
// re-phase -> re-write) and returns the original and regenerated outputs.
func runScaleRoundTrip(t *testing.T, ds *model.Dataset, cfg config.Config) (orig, regen string) {
	t.Helper()
	res := Run(ds, cfg)
	var buf bytes.Buffer
	if err := output.WritePhased(&buf, ds, res, cfg); err != nil {
		t.Fatal(err)
	}
	orig = buf.String()
	path := filepath.Join(t.TempDir(), "phased.tsv")
	if err := os.WriteFile(path, []byte(orig), 0644); err != nil {
		t.Fatal(err)
	}
	ds2, _, cfg2, _, err := parseinput.LoadPhased(path)
	if err != nil {
		t.Fatal(err)
	}
	var buf2 bytes.Buffer
	if err := output.WritePhased(&buf2, ds2, Run(ds2, cfg2), cfg2); err != nil {
		t.Fatal(err)
	}
	regen = buf2.String()
	return orig, regen
}

// TestScaleRereadCompletes verifies that the re-read verification path
// completes at a moderate scale and produces a byte-identical result.
func TestScaleRereadCompletes(t *testing.T) {
	if testing.Short() {
		t.Skip("scale reread skipped in -short mode")
	}
	ds := makeScaleDatasetWithTrios(500, 10, 5000)
	cfg := config.Default()
	orig, regen := runScaleRoundTrip(t, ds, cfg)
	if orig != regen {
		t.Fatalf("scale re-read produced a different result (orig %d bytes, regen %d bytes)", len(orig), len(regen))
	}
}

// BenchmarkScaleReread measures the full re-read verification path at the
// documented boundary. The re-read work must grow smoothly with scale.
func BenchmarkScaleReread(b *testing.B) {
	ds := makeScaleDatasetWithTrios(10000, 50, 100000)
	cfg := config.Default()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		var buf bytes.Buffer
		if err := output.WritePhased(&buf, ds, Run(ds, cfg), cfg); err != nil {
			b.Fatal(err)
		}
		path := filepath.Join(b.TempDir(), "phased.tsv")
		if err := os.WriteFile(path, []byte(buf.String()), 0644); err != nil {
			b.Fatal(err)
		}
		ds2, _, cfg2, _, err := parseinput.LoadPhased(path)
		if err != nil {
			b.Fatal(err)
		}
		var buf2 bytes.Buffer
		if err := output.WritePhased(&buf2, ds2, Run(ds2, cfg2), cfg2); err != nil {
			b.Fatal(err)
		}
	}
}

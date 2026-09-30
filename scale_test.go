package phasing

import (
	"runtime"
	"strconv"
	"strings"
	"testing"
	"time"
)

// Declared supported scale (see docs/ALGORITHM.md):
// 10,000 variant sites, 100 samples, 2,000,000 read observations.
const (
	scaleVariants     = 10_000
	scaleSamples      = 100
	scaleObservations = 2_000_000
	scaleTimeLimit    = 30 * time.Second
)

func TestScaleBoundary(t *testing.T) {
	if testing.Short() {
		t.Skip("scale boundary test skipped in -short mode")
	}
	dir := t.TempDir()

	var vb strings.Builder
	vb.WriteString("variant_id\tchrom\tpos\tref\talt\n")
	for i := 0; i < scaleVariants; i++ {
		vb.WriteString("v" + strconv.Itoa(i) + "\tchr1\t" +
			strconv.Itoa(100+i*100) + "\tA\tT\n")
	}
	vp := writeFile(t, dir, "variants.tsv", vb.String())

	// Half the samples are heterozygous everywhere; reads are generated only
	// for those samples. Dense reads: each read spans 2 adjacent hets, so the
	// target number of observations is achieved with observations/2 reads.
	var gb strings.Builder
	gb.WriteString("sample_id")
	for i := 0; i < scaleVariants; i++ {
		gb.WriteString("\tv" + strconv.Itoa(i))
	}
	gb.WriteString("\n")
	hetSamples := scaleSamples / 2
	for s := 0; s < scaleSamples; s++ {
		gb.WriteString("s" + strconv.Itoa(s))
		for i := 0; i < scaleVariants; i++ {
			if s < hetSamples {
				gb.WriteString("\t0/1")
			} else {
				gb.WriteString("\t0/0")
			}
		}
		gb.WriteString("\n")
	}
	gp := writeFile(t, dir, "genotypes.tsv", gb.String())

	readsPerSample := scaleObservations / 2 / hetSamples
	var rb strings.Builder
	rb.Grow(scaleObservations * 12)
	rb.WriteString("read_id\tsample_id\tvariant_id\tallele\n")
	readID := 0
	for s := 0; s < hetSamples; s++ {
		sample := "s" + strconv.Itoa(s)
		// Pair adjacent sites; alternate phase deterministically.
		for e := 0; e < readsPerSample && 2*e+1 < scaleVariants; e++ {
			v1, v2 := 2*e, 2*e+1
			a1, a2 := int8(0), int8(e%2)
			rb.WriteString("r" + strconv.Itoa(readID) + "\t" + sample +
				"\tv" + strconv.Itoa(v1) + "\t" + strconv.Itoa(int(a1)) + "\n")
			rb.WriteString("r" + strconv.Itoa(readID) + "\t" + sample +
				"\tv" + strconv.Itoa(v2) + "\t" + strconv.Itoa(int(a2)) + "\n")
			readID++
		}
	}
	rp := writeFile(t, dir, "reads.tsv", rb.String())

	start := time.Now()
	in, err := ParseInput(vp, gp, rp, "", "")
	if err != nil {
		t.Fatalf("ParseInput at scale: %v", err)
	}
	res, err := Run(in)
	if err != nil {
		t.Fatalf("Run at scale: %v", err)
	}
	if err := WriteResult(dir+"/out.json", res); err != nil {
		t.Fatalf("WriteResult at scale: %v", err)
	}
	elapsed := time.Since(start)

	if elapsed > scaleTimeLimit {
		t.Fatalf("end-to-end at declared scale took %v, limit %v", elapsed, scaleTimeLimit)
	}
	if res.Stats.VariantSites != scaleVariants || res.Stats.Samples != scaleSamples {
		t.Fatalf("scale stats mismatch: %+v", res.Stats)
	}
	var phased int
	for _, sr := range res.Samples[:hetSamples] {
		phased += countPhased(sr)
	}
	if phased == 0 {
		t.Fatal("expected phased hets at scale")
	}

	// Measure live heap after a forced GC; transient serialization buffers must
	// not be confused with the service's working set. The documented ceiling is
	// on retained memory at the declared boundary.
	runtime.GC()
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	ceiling := uint64(512) * 1024 * 1024
	if m.HeapAlloc > ceiling {
		t.Fatalf("retained heap %d bytes exceeds 512 MiB ceiling", m.HeapAlloc)
	}
	t.Logf("scale OK: elapsed=%v phased_hets=%d heapMiB=%.1f totalAllocMiB=%.1f",
		elapsed, phased, float64(m.HeapAlloc)/(1024*1024), float64(m.TotalAlloc)/(1024*1024))
}

func countPhased(sr SampleResult) int {
	n := 0
	for _, s := range sr.Sites {
		if s.Status == SitePhased && s.H1 != s.H2 {
			n++
		}
	}
	return n
}

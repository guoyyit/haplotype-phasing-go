package parseinput

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/output"
	"github.com/guoyyit/haplotype-phasing-go/internal/phase"
)

func writeTmp(t *testing.T, name, content string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

const variantsTSV = `#variants
v1	chr1	100	A	G
v2	chr1	200	C	T
`

const genotypesTSV = `#genotypes
s1	v1	A/G
s1	v2	C/T
`

const readsTSV = `#reads
r1	s1	v1	A
r1	s1	v2	T
r2	s1	v1	G
r2	s1	v2	C
`

// LoadDataset parses all four inputs and tolerates absent reads/trios.
func TestLoadDatasetBasic(t *testing.T) {
	v := writeTmp(t, "v.tsv", variantsTSV)
	g := writeTmp(t, "g.tsv", genotypesTSV)
	r := writeTmp(t, "r.tsv", readsTSV)
	ds, notices, err := LoadDataset(v, g, r, "")
	if err != nil {
		t.Fatal(err)
	}
	if len(ds.Variants) != 2 || len(ds.Samples) != 1 {
		t.Fatalf("unexpected dataset dims: %d variants %d samples", len(ds.Variants), len(ds.Samples))
	}
	if len(ds.Reads) != 2 || len(ds.Reads[0].Observations) != 2 {
		t.Fatalf("reads not grouped correctly: %+v", ds.Reads)
	}
	_ = notices
}

// Missing/partial genotypes parse without error.
func TestLoadMissingGenotypes(t *testing.T) {
	v := writeTmp(t, "v.tsv", variantsTSV)
	g := writeTmp(t, "g.tsv", "#genotypes\ns1\tv1\t./.\ns1\tv2\tC/.\n")
	ds, _, err := LoadDataset(v, g, "", "")
	if err != nil {
		t.Fatal(err)
	}
	if ds.GT["s1"][0].Ploidy != 0 || ds.GT["s1"][1].Ploidy != 1 {
		t.Fatalf("missing/partial not parsed: %+v", ds.GT["s1"])
	}
}

// Writing then reading a phased file must reproduce identical output.
func TestPhasedRoundTrip(t *testing.T) {
	v := writeTmp(t, "v.tsv", variantsTSV)
	g := writeTmp(t, "g.tsv", genotypesTSV)
	r := writeTmp(t, "r.tsv", readsTSV)
	ds, _, err := LoadDataset(v, g, r, "")
	if err != nil {
		t.Fatal(err)
	}
	cfg := config.Default()
	res := phase.Run(ds, cfg)
	var buf bytes.Buffer
	if err := output.WritePhased(&buf, ds, res, cfg); err != nil {
		t.Fatal(err)
	}
	phasedPath := writeTmp(t, "phased.tsv", buf.String())

	ds2, _, cfg2, _, version, err := LoadPhased(phasedPath)
	if err != nil {
		t.Fatal(err)
	}
	if version != "1.1" {
		t.Fatalf("expected current format version 1.1, got %q", version)
	}
	var buf2 bytes.Buffer
	if err := output.WritePhasedVersion(&buf2, ds2, phase.Run(ds2, cfg2), cfg2, version); err != nil {
		t.Fatal(err)
	}
	if buf.String() != buf2.String() {
		t.Fatal("phased output is not stable across write/re-read/re-write")
	}
}

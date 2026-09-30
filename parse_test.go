package phasing

import (
	"os"
	"path/filepath"
	"testing"
)

func writeFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

const testVariants = `# comment
variant_id	chrom	pos	ref	alt
v2	chr1	200	A	T
v1	chr1	100	C	G
v3	chr2	50	G	A
`

const testGenotypes = `sample_id	v1	v2	v3
child	0/1	1/1	./.
father	0|0	0/1	0/1
mother	1/1	./1	.
unrelated	1/0	./.	0/.
`

const testReads = `read_id	sample_id	variant_id	allele	weight
r1	child	v1	0
r1	child	v2	1
r2	child	v1	1	0.5
`

func TestParseInputOrderingAndMissing(t *testing.T) {
	dir := t.TempDir()
	vp := writeFile(t, dir, "variants.tsv", testVariants)
	gp := writeFile(t, dir, "genotypes.tsv", testGenotypes)
	rp := writeFile(t, dir, "reads.tsv", testReads)

	in, err := ParseInput(vp, gp, rp, "", "")
	if err != nil {
		t.Fatalf("ParseInput: %v", err)
	}
	// Canonical order: chr1:100 (v1), chr1:200 (v2), chr2:50 (v3).
	wantOrder := []string{"v1", "v2", "v3"}
	for i, want := range wantOrder {
		if in.Variants[i].ID != want {
			t.Fatalf("variant[%d] = %s, want %s", i, in.Variants[i].ID, want)
		}
	}

	child := in.Samples[0]
	if child.ID != "child" {
		t.Fatalf("first sample = %q", child.ID)
	}
	// child: v1 0/1 het, v2 1/1 hom, v3 ./. missing.
	if !isHet(child.GT[0]) || !isHomAlt(child.GT[1]) || !isMissing(child.GT[2]) {
		t.Fatalf("unexpected child genotypes: %+v", child.GT)
	}
	// mother v2 ./1 partial deletion.
	mother := in.Samples[2]
	if mother.GT[1].A1 != -1 || mother.GT[1].A2 != 1 {
		t.Fatalf("partial genotype not parsed: %+v", mother.GT[1])
	}

	if len(in.Reads) != 2 || len(in.Reads[0].Obs) != 2 {
		t.Fatalf("reads not grouped correctly: %+v", in.Reads)
	}
	if in.Reads[0].Obs[0].Variant != 0 || in.Reads[0].Obs[1].Variant != 1 {
		t.Fatalf("read obs not in canonical variant order: %+v", in.Reads[0].Obs)
	}
	if in.Reads[1].Obs[0].Weight != 0.5 {
		t.Fatalf("weight not parsed: %+v", in.Reads[1].Obs[0])
	}

	// Defaults resolved.
	if in.Config.MinEdgeSupport != 1 || in.Config.ReadErrorRate != 0.01 {
		t.Fatalf("defaults not applied: %+v", in.Config)
	}
}

func TestParseTrioValidation(t *testing.T) {
	dir := t.TempDir()
	vp := writeFile(t, dir, "variants.tsv", testVariants)
	gp := writeFile(t, dir, "genotypes.tsv", testGenotypes)
	good := writeFile(t, dir, "trio.tsv", "child\tfather\tmother\nchild\tfather\tmother\n")
	in, err := ParseInput(vp, gp, "", good, "")
	if err != nil {
		t.Fatalf("valid trio rejected: %v", err)
	}
	if in.Trio.Child != "child" || in.Trio.Father != "father" || in.Trio.Mother != "mother" {
		t.Fatalf("trio not parsed: %+v", in.Trio)
	}

	bad := writeFile(t, dir, "trio_bad.tsv", "child\tfather\tmother\nchild\tghost\tmother\n")
	if _, err := ParseInput(vp, gp, "", bad, ""); err == nil {
		t.Fatal("unknown father sample should be rejected")
	}
}

func TestParseReadsErrors(t *testing.T) {
	dir := t.TempDir()
	vp := writeFile(t, dir, "variants.tsv", testVariants)
	gp := writeFile(t, dir, "genotypes.tsv", testGenotypes)

	conflict := writeFile(t, dir, "reads.tsv",
		"read_id\tsample_id\tvariant_id\tallele\nr1\tchild\tv1\t0\nr1\tchild\tv1\t1\n")
	if _, err := ParseInput(vp, gp, conflict, "", ""); err == nil {
		t.Fatal("intra-read allele conflict should be rejected")
	}

	multiSample := writeFile(t, dir, "reads2.tsv",
		"read_id\tsample_id\tvariant_id\tallele\nr1\tchild\tv1\t0\nr1\tfather\tv1\t1\n")
	if _, err := ParseInput(vp, gp, multiSample, "", ""); err == nil {
		t.Fatal("read across multiple samples should be rejected")
	}
}

func TestConfigFileOverride(t *testing.T) {
	dir := t.TempDir()
	vp := writeFile(t, dir, "variants.tsv", testVariants)
	gp := writeFile(t, dir, "genotypes.tsv", testGenotypes)
	cp := writeFile(t, dir, "config.json",
		`{"min_edge_support":3,"read_error_rate":0.02,"seed_enabled":true,"random_seed":42}`)
	in, err := ParseInput(vp, gp, "", "", cp)
	if err != nil {
		t.Fatalf("config parse: %v", err)
	}
	if in.Config.MinEdgeSupport != 3 || in.Config.ReadErrorRate != 0.02 ||
		!in.Config.SeedEnabled || in.Config.RandomSeed != 42 {
		t.Fatalf("config overrides not applied: %+v", in.Config)
	}
}

func isHet(c GenotypeCall) bool {
	return c.Present && c.A1 >= 0 && c.A2 >= 0 && c.A1 != c.A2
}

func isHomAlt(c GenotypeCall) bool {
	return c.Present && c.A1 == 1 && c.A2 == 1
}

func isMissing(c GenotypeCall) bool {
	return !c.Present || (c.A1 == -1 && c.A2 == -1)
}

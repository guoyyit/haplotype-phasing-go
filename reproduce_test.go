package phasing

import (
	"bytes"
	"testing"
)

func marshalResult(t *testing.T, r *Result) []byte {
	t.Helper()
	dir := t.TempDir()
	p := dir + "/out.json"
	if err := WriteResult(p, r); err != nil {
		t.Fatalf("WriteResult: %v", err)
	}
	b, err := ReadResultFile(p)
	if err != nil {
		t.Fatalf("read back: %v", err)
	}
	return b
}

// Two runs over identical input/config must emit identical bytes.
func TestDeterministicRepeatedRuns(t *testing.T) {
	variants := `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`
	genotypes := "sample_id\tv1\tv2\ns1\t0/1\t0/1\n"
	// Exactly balanced read evidence forces tie resolution.
	reads := `read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
r2	s1	v1	0
r2	s1	v2	0
`
	run := func() []byte {
		return marshalResult(t, runFixture(t, variants, genotypes, reads))
	}
	first := run()
	for i := 0; i < 3; i++ {
		if !bytes.Equal(first, run()) {
			t.Fatalf("run %d differs from first run", i+1)
		}
	}
}

// Same seed must reproduce; a different seed may resolve ties differently.
func TestSeededTieReproducible(t *testing.T) {
	dir := t.TempDir()
	variants := `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`
	genotypes := "sample_id\tv1\tv2\ns1\t0/1\t0/1\n"
	reads := `read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
r2	s1	v1	0
r2	s1	v2	0
`
	vp := writeFile(t, dir, "variants.tsv", variants)
	gp := writeFile(t, dir, "genotypes.tsv", genotypes)
	rp := writeFile(t, dir, "reads.tsv", reads)

	runWith := func(seed int64) []byte {
		cp := writeFile(t, dir, "config.json",
			`{"seed_enabled":true,"random_seed":`+itoa(seed)+
				`,"min_block_confidence":0.5}`)
		in, err := ParseInput(vp, gp, rp, "", cp)
		if err != nil {
			t.Fatal(err)
		}
		res, err := Run(in)
		if err != nil {
			t.Fatal(err)
		}
		return marshalResult(t, res)
	}

	a := runWith(7)
	b := runWith(7)
	if !bytes.Equal(a, b) {
		t.Fatal("same seed produced different outputs")
	}

	// Across many distinct seeds the PRNG-driven tie must be able to take both
	// values; otherwise the seed setting would have no observable effect.
	sawParity := map[int8]bool{}
	for seed := int64(0); seed < 64; seed++ {
		dir2 := t.TempDir()
		cp2 := writeFile(t, dir2, "config.json",
			`{"seed_enabled":true,"random_seed":`+itoa(seed)+
				`,"min_block_confidence":0.5}`)
		vp2 := writeFile(t, dir2, "variants.tsv", variants)
		gp2 := writeFile(t, dir2, "genotypes.tsv", genotypes)
		rp2 := writeFile(t, dir2, "reads.tsv", reads)
		in2, err := ParseInput(vp2, gp2, rp2, "", cp2)
		if err != nil {
			t.Fatal(err)
		}
		res2, err := Run(in2)
		if err != nil {
			t.Fatal(err)
		}
		sawParity[res2.Samples[0].Sites[1].H1] = true
	}
	if len(sawParity) < 2 {
		t.Fatalf("seeded ties never varied: %v", sawParity)
	}
}

// Rephasing a written result yields byte-identical output.
func TestRephaseRoundTrip(t *testing.T) {
	res := runFixture(t, threeHetVariants,
		"sample_id\tv1\tv2\tv3\ns1\t0/1\t0/1\t0/1\n",
		`read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
r2	s1	v2	1
r2	s1	v3	0
`)
	dir := t.TempDir()
	p := dir + "/out.json"
	if err := WriteResult(p, res); err != nil {
		t.Fatal(err)
	}
	prior, err := ReadResult(p)
	if err != nil {
		t.Fatal(err)
	}

	// Rebuild the same input independently.
	vp := writeFile(t, dir, "variants.tsv", threeHetVariants)
	gp := writeFile(t, dir, "genotypes.tsv", "sample_id\tv1\tv2\tv3\ns1\t0/1\t0/1\t0/1\n")
	rp := writeFile(t, dir, "reads.tsv",
		`read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
r2	s1	v2	1
r2	s1	v3	0
`)
	in, err := ParseInput(vp, gp, rp, "", "")
	if err != nil {
		t.Fatal(err)
	}
	rephased, err := Rephase(in, prior)
	if err != nil {
		t.Fatal(err)
	}
	if err := WriteResult(dir+"/out2.json", rephased); err != nil {
		t.Fatal(err)
	}
	orig, _ := ReadResultFile(p)
	again, _ := ReadResultFile(dir + "/out2.json")
	if !bytes.Equal(orig, again) {
		t.Fatal("rephased output is not byte identical")
	}
}

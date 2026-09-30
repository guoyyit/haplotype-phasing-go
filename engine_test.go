package phasing

import "testing"

func runFixture(t *testing.T, variants, genotypes, reads string) *Result {
	t.Helper()
	dir := t.TempDir()
	vp := writeFile(t, dir, "variants.tsv", variants)
	gp := writeFile(t, dir, "genotypes.tsv", genotypes)
	rp := writeFile(t, dir, "reads.tsv", reads)
	in, err := ParseInput(vp, gp, rp, "", "")
	if err != nil {
		t.Fatalf("ParseInput: %v", err)
	}
	res, err := Run(in)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	return res
}

const threeHetVariants = `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
v3	chr1	300	G	A
`

// One sample, three heterozygous sites; reads consistently say
// v1=0,v2=1 (flip) and v2=1,v3=0 (flip), producing one block with
// h1 alleles 0,1,0 for the anchored orientation.
func TestRunBasicReadPhasing(t *testing.T) {
	genotypes := "sample_id\tv1\tv2\tv3\ns1\t0/1\t0/1\t0/1\n"
	reads := `read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
r2	s1	v2	1
r2	s1	v3	0
`
	res := runFixture(t, threeHetVariants, genotypes, reads)
	if len(res.Samples) != 1 {
		t.Fatalf("samples = %d", len(res.Samples))
	}
	sr := res.Samples[0]
	if len(sr.Blocks) != 1 {
		t.Fatalf("blocks = %d, want 1", len(sr.Blocks))
	}
	wantH1 := []int8{0, 1, 0}
	for v := 0; v < 3; v++ {
		site := sr.Sites[v]
		if site.Status != SitePhased {
			t.Fatalf("site %d status = %s, want phased", v, site.Status)
		}
		if site.H1 != wantH1[v] || site.H2 != 1-wantH1[v] {
			t.Fatalf("site %d = %d/%d, want %d/%d", v, site.H1, site.H2, wantH1[v], 1-wantH1[v])
		}
	}
	b := sr.Blocks[0]
	if b.Hets != 3 || b.ExpectedErrorRate < 0 || b.ExpectedErrorRate > 1 {
		t.Fatalf("bad block info: %+v", b)
	}
	if b.Confidence <= 0 || b.Confidence > 1 {
		t.Fatalf("bad confidence: %v", b.Confidence)
	}
	if res.Stats.PhasedHets != 3 {
		t.Fatalf("phased hets = %d", res.Stats.PhasedHets)
	}
}

// Conflicting reads: two say flip, one says same at v1-v2. The majority wins
// and the losing read must be counted and reported.
func TestRunReadConflictCounted(t *testing.T) {
	genotypes := "sample_id\tv1\tv2\ns1\t0/1\t0/1\n"
	reads := `read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
r2	s1	v1	0
r2	s1	v2	1
r3	s1	v1	0
r3	s1	v2	0
`
	res := runFixture(t,
		`variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`, genotypes, reads)
	if res.Conflicts.ReadReadConflicts != 1 {
		t.Fatalf("read_read conflicts = %d, want 1", res.Conflicts.ReadReadConflicts)
	}
	if len(res.Conflicts.Entries) != 1 || res.Conflicts.Entries[0].ReadID != "r3" {
		t.Fatalf("conflict entry wrong: %+v", res.Conflicts.Entries)
	}
}

// No evidence at v2-v3 forces a block break: two singleton unphased regions
// split the phasing even though v1-v2 has one supporting read.
func TestRunBlockBreakOnInsufficientEvidence(t *testing.T) {
	genotypes := "sample_id\tv1\tv2\tv3\ns1\t0/1\t0/1\t0/1\n"
	reads := `read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
`
	res := runFixture(t, threeHetVariants, genotypes, reads)
	sr := res.Samples[0]
	if len(sr.Blocks) != 1 {
		t.Fatalf("blocks = %d, want 1 (v1-v2 only)", len(sr.Blocks))
	}
	if sr.Sites[2].Status != SiteUnphased {
		t.Fatalf("v3 status = %s, want unphased", sr.Sites[2].Status)
	}
}

// High support threshold rejects a single-read edge and breaks the block.
func TestRunMinEdgeSupportThreshold(t *testing.T) {
	dir := t.TempDir()
	vp := writeFile(t, dir, "variants.tsv",
		`variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`)
	gp := writeFile(t, dir, "genotypes.tsv", "sample_id\tv1\tv2\ns1\t0/1\t0/1\n")
	rp := writeFile(t, dir, "reads.tsv",
		`read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
`)
	cp := writeFile(t, dir, "config.json", `{"min_edge_support":2}`)
	in, err := ParseInput(vp, gp, rp, "", cp)
	if err != nil {
		t.Fatalf("ParseInput: %v", err)
	}
	res, err := Run(in)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	sr := res.Samples[0]
	if len(sr.Blocks) != 0 {
		t.Fatalf("blocks = %d, want 0 under min_edge_support=2", len(sr.Blocks))
	}
	if sr.Sites[0].Status != SiteUnphased || sr.Sites[1].Status != SiteUnphased {
		t.Fatalf("sites should be unphased: %+v", sr.Sites)
	}
}

// A read observing an alt allele at a 0/0 site contradicts the genotype and
// must be counted rather than silently dropped or used as evidence.
func TestRunReadGenotypeConflictCounted(t *testing.T) {
	genotypes := "sample_id\tv1\tv2\ns1\t0/0\t0/1\n"
	reads := `read_id	sample_id	variant_id	allele
r1	s1	v1	1
r1	s1	v2	1
`
	res := runFixture(t,
		`variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`, genotypes, reads)
	if res.Conflicts.ReadGenotypeConflicts != 1 {
		t.Fatalf("read_genotype conflicts = %d, want 1",
			res.Conflicts.ReadGenotypeConflicts)
	}
	if len(res.Conflicts.Entries) != 1 ||
		res.Conflicts.Entries[0].Kind != "read_genotype" {
		t.Fatalf("conflict entry wrong: %+v", res.Conflicts.Entries)
	}
}

// A single unanimous read has posterior 1-e (0.99 by default). Requiring
// 0.999 confidence rejects it even though the vote count is sufficient, so the
// block must break deterministically.
func TestRunMinConfidenceThreshold(t *testing.T) {
	dir := t.TempDir()
	vp := writeFile(t, dir, "variants.tsv",
		`variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`)
	gp := writeFile(t, dir, "genotypes.tsv", "sample_id\tv1\tv2\ns1\t0/1\t0/1\n")
	rp := writeFile(t, dir, "reads.tsv",
		`read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
`)
	cp := writeFile(t, dir, "config.json", `{"min_block_confidence":0.999}`)
	in, err := ParseInput(vp, gp, rp, "", cp)
	if err != nil {
		t.Fatalf("ParseInput: %v", err)
	}
	res, err := Run(in)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	sr := res.Samples[0]
	if len(sr.Blocks) != 0 {
		t.Fatalf("blocks = %d, want 0 under 0.999 confidence floor", len(sr.Blocks))
	}
}

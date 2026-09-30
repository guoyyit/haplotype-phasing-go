package phasing

import "testing"

func runTrioFixture(t *testing.T, variants, genotypes, reads, pedigree string, config ...string) *Result {
	t.Helper()
	dir := t.TempDir()
	vp := writeFile(t, dir, "variants.tsv", variants)
	gp := writeFile(t, dir, "genotypes.tsv", genotypes)
	var rp string
	if reads != "" {
		rp = writeFile(t, dir, "reads.tsv", reads)
	}
	pp := writeFile(t, dir, "pedigree.tsv", pedigree)
	var cp string
	if len(config) > 0 {
		cp = writeFile(t, dir, "config.json", config[0])
	}
	in, err := ParseInput(vp, gp, rp, pp, cp)
	if err != nil {
		t.Fatalf("ParseInput: %v", err)
	}
	res, err := Run(in)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	return res
}

const trioPedigree = "child\tfather\tmother\nchild\tfather\tmother\n"

// Two het sites: father is 0/0 at both (transmits 0), mother 1/1 at both
// (transmits 1). The child is 0/1 at both and must inherit h1=0 from father
// and h2=1 from mother, with both sites in one pedigree-anchored block even
// without any reads.
func TestTrioInheritanceBlock(t *testing.T) {
	variants := `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`
	genotypes := "sample_id\tv1\tv2\nchild\t0/1\t0/1\nfather\t0/0\t0/0\nmother\t1/1\t1/1\n"
	res := runTrioFixture(t, variants, genotypes, "", trioPedigree)
	sr := sampleResult(t, res, "child")
	if len(sr.Blocks) != 1 {
		t.Fatalf("blocks = %d, want 1 trio block", len(sr.Blocks))
	}
	for v := 0; v < 2; v++ {
		s := sr.Sites[v]
		if !s.OriginResolved {
			t.Fatalf("site %d origin not resolved", v)
		}
		if !s.PaternalH1 {
			t.Fatalf("site %d h1 not marked paternal", v)
		}
		if s.H1 != 0 || s.H2 != 1 {
			t.Fatalf("site %d = %d/%d, want 0/1 paternal/maternal", v, s.H1, s.H2)
		}
	}
	if res.Conflicts.TrioMendelViolations != 0 {
		t.Fatalf("mendel violations = %d, want 0", res.Conflicts.TrioMendelViolations)
	}
}

// Child 1/1 with father 0/0 is impossible and must be flagged and excluded
// from correct phasing counts.
func TestTrioMendelViolationFlagged(t *testing.T) {
	variants := `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`
	// v1 consistent (child 0/1, father 0/0, mother 1/1); v2 impossible.
	genotypes := "sample_id\tv1\tv2\nchild\t0/1\t1/1\nfather\t0/0\t0/0\nmother\t1/1\t0/1\n"
	res := runTrioFixture(t, variants, genotypes, "", trioPedigree)
	sr := sampleResult(t, res, "child")
	if !sr.Sites[1].MendelViolation {
		t.Fatalf("site v2 should be flagged Mendel violation: %+v", sr.Sites[1])
	}
	if sr.Sites[1].Status != SiteUnphased {
		t.Fatalf("violating site must be unphased, got %s", sr.Sites[1].Status)
	}
	if res.Conflicts.TrioMendelViolations != 1 {
		t.Fatalf("violations = %d, want 1", res.Conflicts.TrioMendelViolations)
	}
	// Only the consistent het counts as correctly phased.
	if res.Stats.CorrectlyPhased < 1 {
		t.Fatalf("correctly phased = %d, want >= 1", res.Stats.CorrectlyPhased)
	}
}

// Read evidence contradicts pedigree: reads claim flip while father/mother
// transmission requires same. Pedigree wins and the conflict is counted.
func TestTrioReadConflictAdjudication(t *testing.T) {
	variants := `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`
	genotypes := "sample_id\tv1\tv2\nchild\t0/1\t0/1\nfather\t0/0\t0/0\nmother\t1/1\t1/1\n"
	// Reads say v1=0,v2=1: same molecule carries 0 and 1, i.e. alleles
	// differ -> flip parity. Pedigree requires both paternal alleles 0 -> same.
	reads := `read_id	sample_id	variant_id	allele
r1	child	v1	0
r1	child	v2	1
`
	res := runTrioFixture(t, variants, genotypes, reads, trioPedigree)
	if res.Conflicts.ReadTrioConflicts != 1 {
		t.Fatalf("read_trio conflicts = %d, want 1", res.Conflicts.ReadTrioConflicts)
	}
	sr := sampleResult(t, res, "child")
	// Pedigree parity (same) kept: h1 = 0 at both sites.
	if sr.Sites[0].H1 != 0 || sr.Sites[1].H1 != 0 {
		t.Fatalf("trio parity not enforced: %d,%d", sr.Sites[0].H1, sr.Sites[1].H1)
	}
}

// Missing/partial genotypes in parents must not fail the run.
func TestTrioMissingParents(t *testing.T) {
	variants := `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`
	genotypes := "sample_id\tv1\tv2\nchild\t0/1\t0/1\nfather\t./.\t./1\nmother\t1/1\t0/.\n"
	res := runTrioFixture(t, variants, genotypes, "", trioPedigree)
	sr := sampleResult(t, res, "child")
	if sr.Sites[0].MendelViolation {
		t.Fatalf("unknown father at v1 cannot be a violation")
	}
	// v1 mother 1/1 transmits 1; child 0/1 -> paternal transmission is 0.
	if !sr.Sites[0].OriginResolved {
		t.Fatalf("v1 origin should be resolved via mother 1/1")
	}
}

func sampleResult(t *testing.T, res *Result, id string) SampleResult {
	t.Helper()
	for _, s := range res.Samples {
		if s.ID == id {
			return s
		}
	}
	t.Fatalf("sample %q not found", id)
	return SampleResult{}
}

package phasing

import "testing"

// Missing, partially missing and homozygous genotypes must not abort the run;
// only fully-known heterozygous sites are eligible for phasing.
func TestMissingAndPartialGenotypes(t *testing.T) {
	variants := `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
v3	chr1	300	G	A
`
	genotypes := "sample_id\tv1\tv2\tv3\ns1\t./.\t./1\t0/0\n"
	reads := `read_id	sample_id	variant_id	allele
r1	s1	v1	0
r1	s1	v2	1
`
	res := runFixture(t, variants, genotypes, reads)
	sr := res.Samples[0]
	if sr.Sites[0].Status != SiteMissing || sr.Sites[0].H1 != -1 {
		t.Fatalf("fully missing site wrong: %+v", sr.Sites[0])
	}
	// v2 partial: known allele 1 preserved on h1, unknown on h2, unphased.
	if sr.Sites[1].Status != SiteUnphased || sr.Sites[1].H1 != 1 || sr.Sites[1].H2 != -1 {
		t.Fatalf("partial site wrong: %+v", sr.Sites[1])
	}
	// v3 hom ref still resolvable and reported phased (trivial homozygote).
	if sr.Sites[2].Status != SitePhased || sr.Sites[2].H1 != 0 || sr.Sites[2].H2 != 0 {
		t.Fatalf("hom ref site wrong: %+v", sr.Sites[2])
	}
	if len(sr.Blocks) != 0 {
		t.Fatalf("missing/partial sites cannot form a block: %+v", sr.Blocks)
	}
	if res.Stats.MissingSites != 1 {
		t.Fatalf("missing sites = %d, want 1", res.Stats.MissingSites)
	}
	if res.Stats.PartialSites != 1 {
		t.Fatalf("partial sites = %d, want 1", res.Stats.PartialSites)
	}
}

// A site with coverage below min_reads_per_site is flagged but must remain
// readable and phasable through pedigree evidence.
func TestLowCoverageFlaggedButProcessed(t *testing.T) {
	dir := t.TempDir()
	variants := `variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
`
	genotypes := "sample_id\tv1\tv2\nchild\t0/1\t0/1\nfather\t0/0\t0/0\nmother\t1/1\t1/1\n"
	ped := "child\tfather\tmother\nchild\tfather\tmother\n"
	cfg := `{"min_reads_per_site":3}`
	vp := writeFile(t, dir, "variants.tsv", variants)
	gp := writeFile(t, dir, "genotypes.tsv", genotypes)
	pp := writeFile(t, dir, "pedigree.tsv", ped)
	cp := writeFile(t, dir, "config.json", cfg)
	in, err := ParseInput(vp, gp, "", pp, cp)
	if err != nil {
		t.Fatalf("ParseInput: %v", err)
	}
	res, err := Run(in)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	child := sampleResult(t, res, "child")
	for v := 0; v < 2; v++ {
		if !child.Sites[v].LowCoverage {
			t.Fatalf("child site %d should be flagged low coverage", v)
		}
		if child.Sites[v].Status != SitePhased {
			t.Fatalf("low coverage trio-resolved site %d should still be phased", v)
		}
	}
	if res.Stats.LowCoverageSites == 0 {
		t.Fatal("low coverage sites must be counted")
	}
	_ = res
}

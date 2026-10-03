package parseinput

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
	"github.com/guoyyit/haplotype-phasing-go/internal/output"
	"github.com/guoyyit/haplotype-phasing-go/internal/phase"
)

// The missing marker is handled on either allele position, and invalid
// alleles are rejected explicitly instead of becoming third alleles.
func TestParseGenotypePositionsAndRejection(t *testing.T) {
	v := model.Variant{ID: "v1", Ref: "A", Alt: "G"}
	cases := []struct {
		token   string
		ploidy  int
		a1      string
		wantErr bool
	}{
		{"./.", 0, "", false},
		{"A/.", 1, "A", false},
		{"./A", 1, "A", false},
		{"G|.", 1, "G", false},
		{".|G", 1, "G", false},
		{"A/G", 2, "A", false},
		{"G/A", 2, "G", false},
		{"A/X", 0, "", true},
		{"X/A", 0, "", true},
		{"./X", 0, "", true},
		{"X/.", 0, "", true},
	}
	for _, tc := range cases {
		gt, err := parseGenotype(tc.token, v)
		if tc.wantErr {
			if err == nil {
				t.Fatalf("genotype %q must be rejected", tc.token)
			}
			if !strings.Contains(err.Error(), "neither ref") {
				t.Fatalf("genotype %q rejection must explain the allele problem: %v", tc.token, err)
			}
			continue
		}
		if err != nil {
			t.Fatalf("genotype %q unexpectedly rejected: %v", tc.token, err)
		}
		if gt.Ploidy != tc.ploidy {
			t.Fatalf("genotype %q: want ploidy %d, got %d", tc.token, tc.ploidy, gt.Ploidy)
		}
		if tc.ploidy >= 1 && gt.A1 != tc.a1 {
			t.Fatalf("genotype %q: want observed allele %q, got %q", tc.token, tc.a1, gt.A1)
		}
	}
}

// End-to-end loading rejects non ref/alt genotype alleles with a locatable
// error naming file, line and allele.
func TestLoadDatasetRejectsThirdAllele(t *testing.T) {
	v := writeTmp(t, "v.tsv", variantsTSV)
	for _, gt := range []string{"s1\tv1\tA/X\n", "s1\tv1\tX/A\n", "s1\tv1\t./X\n"} {
		g := writeTmp(t, "g.tsv", "#genotypes\n"+gt)
		_, _, err := LoadDataset(v, g, "", "")
		if err == nil || !strings.Contains(err.Error(), "neither ref") {
			t.Fatalf("genotype %q must fail with an explicit allele error, got %v", gt, err)
		}
	}
}

// Partial genotypes with the missing marker on either side parse to
// ploidy-1 observations and survive the full round trip.
func TestPartialMissingRoundTripBothSides(t *testing.T) {
	v := writeTmp(t, "v.tsv", variantsTSV)
	g := writeTmp(t, "g.tsv", "#genotypes\ns1\tv1\tA/.\ns1\tv2\t./T\n")
	ds, _, err := LoadDataset(v, g, "", "")
	if err != nil {
		t.Fatal(err)
	}
	if ds.GT["s1"][0].Ploidy != 1 || ds.GT["s1"][0].A1 != "A" {
		t.Fatalf("right-side partial parsed wrong: %+v", ds.GT["s1"][0])
	}
	if ds.GT["s1"][1].Ploidy != 1 || ds.GT["s1"][1].A1 != "T" {
		t.Fatalf("left-side partial parsed wrong: %+v", ds.GT["s1"][1])
	}
	cfg := config.Default()
	res := phase.Run(ds, cfg)
	if res.Samples["s1"].Sites[0].Phased || res.Samples["s1"].Sites[1].Phased {
		t.Fatal("partial sites must never be PHASED")
	}
	assertRoundTrip(t, ds, cfg)
}

// Files produced by the legacy 1.0 writer, including rows for violating
// heterozygous sites (encoded as two missing calls), must still verify
// byte-for-byte under the current service.
func TestLegacyV1RereadCompat(t *testing.T) {
	cfg := config.Default()
	ds := &model.Dataset{
		Variants: []model.Variant{{ID: "v1", Ref: "A", Alt: "G"}, {ID: "v2", Ref: "C", Alt: "T"}},
		Samples:  []string{"c", "f", "m"},
		GT: map[string][]model.Genotype{
			"c": {{Present: true, Ploidy: 2, A1: "A", A2: "G"}, {Present: true, Ploidy: 2, A1: "T", A2: "T"}},
			"f": {{Present: true, Ploidy: 2, A1: "A", A2: "A"}, {Present: true, Ploidy: 2, A1: "C", A2: "C"}},
			"m": {{Present: true, Ploidy: 2, A1: "A", A2: "A"}, {Present: true, Ploidy: 2, A1: "C", A2: "C"}},
		},
		Trios: []model.Trio{{Child: "c", Father: "f", Mother: "m"}},
	}
	res := phase.Run(ds, cfg)
	var legacy strings.Builder
	if err := output.WritePhasedVersion(&legacy, ds, res, cfg, "1.0"); err != nil {
		t.Fatal(err)
	}
	body := legacy.String()
	if !strings.Contains(body, "#haplophased\t1.0\n") {
		t.Fatal("legacy writer must emit the 1.0 version header")
	}
	// Both violation rows exist; the het one uses the 1.0 dotted encoding.
	if !strings.Contains(body, "c\t0\tv1\tUNPHASED\t1\tA\tG\t.\t.\t1") {
		t.Fatal("1.0 violating het row must keep the legacy dotted encoding")
	}
	path := filepath.Join(t.TempDir(), "legacy.tsv")
	if err := os.WriteFile(path, []byte(body), 0644); err != nil {
		t.Fatal(err)
	}
	ds2, _, cfg2, _, version, err := LoadPhased(path)
	if err != nil {
		t.Fatal(err)
	}
	if version != "1.0" {
		t.Fatalf("expected version 1.0, got %q", version)
	}
	var regen strings.Builder
	if err := output.WritePhasedVersion(&regen, ds2, phase.Run(ds2, cfg2), cfg2, version); err != nil {
		t.Fatal(err)
	}
	if body != regen.String() {
		t.Fatal("legacy 1.0 file must re-emit byte-for-byte")
	}
}

// A format version the service does not know must be rejected with a
// locatable parse error rather than an unexplained mismatch later.
func TestUnknownVersionRejected(t *testing.T) {
	v := writeTmp(t, "v.tsv", variantsTSV)
	g := writeTmp(t, "g.tsv", genotypesTSV)
	ds, _, err := LoadDataset(v, g, "", "")
	if err != nil {
		t.Fatal(err)
	}
	var body strings.Builder
	_ = output.WritePhased(&body, ds, phase.Run(ds, config.Default()), config.Default())
	future := strings.Replace(body.String(), "#haplophased\t1.1", "#haplophased\t9.9", 1)
	path := filepath.Join(t.TempDir(), "future.tsv")
	if err := os.WriteFile(path, []byte(future), 0644); err != nil {
		t.Fatal(err)
	}
	if _, _, _, _, _, err := LoadPhased(path); err == nil ||
		!strings.Contains(err.Error(), "unsupported phased format version") {
		t.Fatalf("unknown version must be rejected explicitly, got %v", err)
	}
}

// A 1.0 file whose row was produced by the historical left-missing-marker
// parsing defect (missing token mistaken for a real allele) must still be
// re-emitted byte-identically: the legacy row content is preserved verbatim
// rather than silently rewritten or rejected.
func TestLegacyV1DefectivePartialRowPreserved(t *testing.T) {
	const file = "#haplophased\t1.0\n" +
		"#variants\n" +
		"variant_id\tchrom\tpos\tref\talt\n" +
		"v1\tc\t1\tA\tG\n" +
		"v2\tc\t2\tC\tT\n" +
		"#trios\n" +
		"child\tfather\tmother\n" +
		"#reads\n" +
		"read_id\tsample\tvariant_id\tallele\n" +
		"#config\n" +
		`{"min_edge_support":2,"min_edge_strength":0.8,"min_block_confidence":0.95,"block_max_error_rate":0.25,"min_site_depth":3,"seed":1,"enable_trio":true,"enable_reads":true}` + "\n" +
		"#phase\n" +
		"sample\tvariant_index\tvariant_id\tstatus\tlow_cov\ta0\ta1\th1\th2\tviolation\n" +
		"s1\t0\tv1\tPHASED\t1\tA\tG\tA0\tA1\t0\n" +
		"s1\t1\tv2\tPHASED\t1\t.\tT\tA0\tA1\t0\n" +
		"#blocks\n" +
		"sample\tblock_id\tstart_idx\tend_idx\tsite_count\thets_count\tsupport\tconflict\tconfidence\terror_rate\tstatus\n" +
		"s1\tBLK-s1-0001\t0\t0\t1\t1\t0\t0\t1.0000\t0.0000\tSINGLETON\n" +
		"s1\tBLK-s1-0002\t1\t1\t1\t1\t0\t0\t1.0000\t0.0000\tSINGLETON\n" +
		"#stats\n" +
		"samples\t1\n" +
		"variants\t2\n" +
		"hets_total\t2\n" +
		"phased_hets\t2\n" +
		"unphased_hets\t0\n" +
		"blocks_total\t2\n" +
		"read_observations\t0\n" +
		"read_pairs\t0\n" +
		"read_edge_supports\t0\n" +
		"read_edge_conflicts\t0\n" +
		"read_conflicts_adjudicated\t0\n" +
		"trio_count\t0\n" +
		"trio_constraint_sites\t0\n" +
		"trio_violations\t0\n" +
		"breaks_by_threshold\t0\n" +
		"low_coverage_sites\t2\n" +
		"missing_genotypes\t0\n"
	path := filepath.Join(t.TempDir(), "legacy-defect.tsv")
	if err := os.WriteFile(path, []byte(file), 0644); err != nil {
		t.Fatal(err)
	}
	ds, _, cfg, _, version, err := LoadPhased(path)
	if err != nil {
		t.Fatal(err)
	}
	var regen strings.Builder
	if err := output.WritePhasedVersion(&regen, ds, phase.Run(ds, cfg), cfg, version); err != nil {
		t.Fatal(err)
	}
	if file != regen.String() {
		t.Fatalf("historical 1.0 defect row must be preserved byte-for-byte")
	}
}

// assertRoundTrip writes the result, reloads it and asserts byte equality.
func assertRoundTrip(t *testing.T, ds *model.Dataset, cfg config.Config) {
	t.Helper()
	res := phase.Run(ds, cfg)
	var first strings.Builder
	if err := output.WritePhased(&first, ds, res, cfg); err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(t.TempDir(), "phased.tsv")
	if err := os.WriteFile(path, []byte(first.String()), 0644); err != nil {
		t.Fatal(err)
	}
	ds2, _, cfg2, _, version, err := LoadPhased(path)
	if err != nil {
		t.Fatal(err)
	}
	var second strings.Builder
	if err := output.WritePhasedVersion(&second, ds2, phase.Run(ds2, cfg2), cfg2, version); err != nil {
		t.Fatal(err)
	}
	if first.String() != second.String() {
		t.Fatal("write/read/re-phase/write is not byte-stable")
	}
}

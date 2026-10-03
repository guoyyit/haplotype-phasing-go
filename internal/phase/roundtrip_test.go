package phase_test

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
	"github.com/guoyyit/haplotype-phasing-go/internal/output"
	"github.com/guoyyit/haplotype-phasing-go/internal/parseinput"
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

// gt2 builds a diploid genotype.
func gt2(a, b string) model.Genotype {
	return model.Genotype{Present: true, Ploidy: 2, A1: a, A2: b}
}

func gt1(a string) model.Genotype {
	return model.Genotype{Present: true, Ploidy: 1, A1: a}
}

func gt0() model.Genotype { return model.Genotype{Present: true, Ploidy: 0} }

func scenarioVariants() []model.Variant {
	return []model.Variant{
		{ID: "v1", Chrom: "c", Pos: 1, Ref: "A", Alt: "G"},
		{ID: "v2", Chrom: "c", Pos: 2, Ref: "C", Alt: "T"},
		{ID: "v3", Chrom: "c", Pos: 3, Ref: "G", Alt: "A"},
		{ID: "v4", Chrom: "c", Pos: 4, Ref: "T", Alt: "C"},
		{ID: "v5", Chrom: "c", Pos: 5, Ref: "A", Alt: "G"},
		{ID: "v6", Chrom: "c", Pos: 6, Ref: "C", Alt: "T"},
	}
}

// scenarioDataset covers every required edge case at once:
//   - child v1 A/G under A/A x A/A parents: heterozygous Mendelian violation
//   - child v2 T/T under C/C x C/C parents: homozygous Mendelian violation
//   - child v3/v4 hets with uniquely forced paternal bits (read/trio edge)
//   - child v5 A/. : missing marker on the right (partial)
//   - child v6 ./T : missing marker on the left (partial)
//   - father v6 is completely missing (. / .)
func scenarioDataset() *model.Dataset {
	v := scenarioVariants()
	child := []model.Genotype{
		gt2("A", "G"), gt2("T", "T"), gt2("G", "A"), gt2("T", "C"), gt1("A"), gt1("T"),
	}
	father := []model.Genotype{
		gt2("A", "A"), gt2("C", "C"), gt2("A", "A"), gt2("T", "T"), gt2("A", "A"), gt0(),
	}
	mother := []model.Genotype{
		gt2("A", "A"), gt2("C", "C"), gt2("G", "G"), gt2("C", "C"), gt2("G", "G"), gt2("C", "C"),
	}
	return &model.Dataset{
		Variants: v,
		Samples:  []string{"child", "father", "mother"},
		GT: map[string][]model.Genotype{
			"child": child, "father": father, "mother": mother,
		},
		Trios: []model.Trio{{Child: "child", Father: "father", Mother: "mother"}},
	}
}

// conflictingReads returns read evidence whose orientation disagrees with
// the trio-forced paternal orientation across the adjacent forced hets
// v3-v4 (paternal v3=A bit1, v4=T bit0 => trio says opp; reads say same),
// plus unanimous supporting reads so the edge connects.
func conflictingReads() []model.Read {
	return []model.Read{
		{ID: "r1", Sample: "child", Observations: []model.ReadObservation{
			{VariantIdx: 2, Allele: "A"}, {VariantIdx: 3, Allele: "C"},
		}},
		{ID: "r2", Sample: "child", Observations: []model.ReadObservation{
			{VariantIdx: 2, Allele: "A"}, {VariantIdx: 3, Allele: "C"},
		}},
		{ID: "r3", Sample: "child", Observations: []model.ReadObservation{
			{VariantIdx: 2, Allele: "A"}, {VariantIdx: 3, Allele: "C"},
		}},
	}
}

// writeRereadRewrite runs the full write -> LoadPhased -> re-phase -> write
// cycle and asserts byte equality plus equal statistics and per-site rows.
func writeRereadRewrite(t *testing.T, ds *model.Dataset, cfg config.Config) {
	t.Helper()
	res := phase.Run(ds, cfg)
	var first bytes.Buffer
	if err := output.WritePhased(&first, ds, res, cfg); err != nil {
		t.Fatal(err)
	}

	// Guard: no row may print the missing marker as a canonical allele.
	for _, line := range bytes.Split(first.Bytes(), []byte("\n")) {
		fields := bytes.Split(line, []byte("\t"))
		if len(fields) == 10 && string(fields[2]) != "variant_id" {
			idx, convErr := strconv.Atoi(string(fields[1]))
			if convErr != nil {
				t.Fatalf("non-numeric variant index in output: %q", fields[1])
			}
			if idx >= 0 && idx < len(ds.Variants) {
				if string(fields[5]) != ds.Variants[idx].Ref ||
					string(fields[6]) != ds.Variants[idx].Alt {
					t.Fatalf("canonical allele columns must equal variant ref/alt: %s", line)
				}
			}
		}
	}

	path := writeTmp(t, "phased.tsv", first.String())
	ds2, res2, cfg2, _, version, err := parseinput.LoadPhased(path)
	if err != nil {
		t.Fatal(err)
	}
	if version != "1.1" {
		t.Fatalf("expected format version 1.1, got %q", version)
	}
	rephased := phase.Run(ds2, cfg2)
	var second bytes.Buffer
	if err := output.WritePhasedVersion(&second, ds2, rephased, cfg2, version); err != nil {
		t.Fatal(err)
	}
	if first.String() != second.String() {
		a, b := splitLines(first.Bytes()), splitLines(second.Bytes())
		for i := 0; i < len(a) || i < len(b); i++ {
			var la, lb []byte
			if i < len(a) {
				la = a[i]
			}
			if i < len(b) {
				lb = b[i]
			}
			if !bytes.Equal(la, lb) {
				t.Fatalf("byte mismatch at line %d:\n original=%q\nrephased=%q", i+1, la, lb)
			}
		}
	}
	if res2.Stats != rephased.Stats {
		t.Fatalf("run stats changed across reread:\noriginal=%+v\nrephased=%+v", res2.Stats, rephased.Stats)
	}
	for _, sample := range res.Order {
		sa, sb := res.Samples[sample], rephased.Samples[sample]
		if sa.HetsTotal != sb.HetsTotal || sa.PhasedHets != sb.PhasedHets ||
			sa.UnphasedHets != sb.UnphasedHets || sa.ViolationHets != sb.ViolationHets {
			t.Fatalf("sample %s counters changed across reread: %+v vs %+v", sample, sa, sb)
		}
		for i := range sa.Sites {
			x, y := sa.Sites[i], sb.Sites[i]
			if x != y {
				t.Fatalf("sample %s site %d changed across reread:\noriginal=%+v\nrephased=%+v", sample, i, x, y)
			}
		}
	}
}

func splitLines(s []byte) [][]byte { return bytes.Split(s, []byte("\n")) }

// All required input combinations must verify byte-for-byte across a full
// write / re-read / re-phase / rewrite cycle.
func TestRereadMatrix(t *testing.T) {
	cases := []struct {
		name      string
		mutate    func(ds *model.Dataset)
		configure func(cfg *config.Config)
	}{
		{"trio_violations", nil, nil},
		{"read_trio_conflict", func(ds *model.Dataset) {
			ds.Reads = conflictingReads()
		}, nil},
		{"partial_and_full_missing", nil, nil},
		{"low_coverage", func(ds *model.Dataset) {
			ds.Reads = conflictingReads()
		}, func(cfg *config.Config) { cfg.MinSiteDepth = 10 }},
		{"no_reads", nil, func(cfg *config.Config) { cfg.EnableReads = false }},
		{"no_trio", func(ds *model.Dataset) { ds.Trios = nil }, nil},
		{"reads_disabled", nil, func(cfg *config.Config) { cfg.EnableReads = false }},
		{"trio_disabled", nil, func(cfg *config.Config) { cfg.EnableTrio = false }},
		{"non_default_config", func(ds *model.Dataset) {
			ds.Reads = conflictingReads()
		}, func(cfg *config.Config) {
			cfg.MinEdgeSupport = 1
			cfg.MinEdgeStrength = 0.5
			cfg.MinBlockConfidence = 0.5
			cfg.BlockMaxErrorRate = 1
			cfg.MinSiteDepth = 2
		}},
		{"non_default_seed", nil, func(cfg *config.Config) { cfg.Seed = 987654321 }},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			ds := scenarioDataset()
			if tc.mutate != nil {
				tc.mutate(ds)
			}
			cfg := config.Default()
			if tc.configure != nil {
				tc.configure(&cfg)
			}
			writeRereadRewrite(t, ds, cfg)
		})
	}
}

// Violating heterozygous and homozygous sites keep their violation flag,
// UNPHASED status, observed alleles and statistics across reread.
func TestViolationRowsRoundTrip(t *testing.T) {
	ds := scenarioDataset()
	cfg := config.Default()
	res := phase.Run(ds, cfg)
	child := res.Samples["child"]

	het := child.Sites[0]
	if !het.Violation || !het.Hetero || het.Phased {
		t.Fatalf("v1 must be a violating unphased heterozygote: %+v", het)
	}
	if het.H1 != model.AlleleA0 || het.H2 != model.AlleleA1 {
		t.Fatalf("violating het must retain observed canonical alleles, got %v/%v", het.H1, het.H2)
	}
	hom := child.Sites[1]
	if !hom.Violation || hom.Hetero || hom.H1 != model.AlleleA1 || hom.H2 != model.AlleleA1 {
		t.Fatalf("v2 must be a violating homozygote carrying A1: %+v", hom)
	}
	if res.Stats.TrioViolations != 2 {
		t.Fatalf("expected 2 trio violations (het + homo), got %d", res.Stats.TrioViolations)
	}
	if child.ViolationHets != 1 {
		t.Fatalf("expected exactly 1 violating het, got %d", child.ViolationHets)
	}
	writeRereadRewrite(t, ds, cfg)

	// Status tokens in the written file must match the documented semantics.
	var buf bytes.Buffer
	_ = output.WritePhased(&buf, ds, res, cfg)
	assertRowStatus(t, buf.Bytes(), "child", 0, model.SiteStatusUnphased, "1")
	assertRowStatus(t, buf.Bytes(), "child", 1, model.SiteStatusUnphased, "1")
	assertRowStatus(t, buf.Bytes(), "child", 4, model.SiteStatusPartial, "0")
	assertRowStatus(t, buf.Bytes(), "child", 5, model.SiteStatusPartial, "0")
	assertRowStatus(t, buf.Bytes(), "father", 5, model.SiteStatusMissing, "0")
}

// Partial genotypes never become PHASED regardless of which side the
// missing marker occupies, and the observed allele is retained.
func TestPartialMissingBothSides(t *testing.T) {
	ds := scenarioDataset()
	res := phase.Run(ds, config.Default())
	child := res.Samples["child"]
	right := child.Sites[4] // v5 A/.
	left := child.Sites[5]  // v6 ./T
	if right.Phased || right.H1 != model.AlleleA0 || right.H2 != model.AlleleMissing {
		t.Fatalf("right-side partial mishandled: %+v", right)
	}
	if left.Phased || left.H1 != model.AlleleA1 || left.H2 != model.AlleleMissing {
		t.Fatalf("left-side partial mishandled: %+v", left)
	}
	fatherMissing := res.Samples["father"].Sites[5]
	if fatherMissing.Phased || fatherMissing.H1 != model.AlleleMissing || fatherMissing.H2 != model.AlleleMissing {
		t.Fatalf("fully missing site mishandled: %+v", fatherMissing)
	}
	if res.Stats.MissingGenotypes < 1 {
		t.Fatal("fully missing genotype must be counted")
	}
	writeRereadRewrite(t, ds, config.Default())
}

// Read evidence that disagrees with the trio-forced orientation is counted
// as a conflict and the whole run still round-trips byte-for-byte.
func TestReadTrioConflictRoundTrip(t *testing.T) {
	ds := scenarioDataset()
	ds.Reads = conflictingReads()
	cfg := config.Default()
	cfg.MinEdgeSupport = 1 // tie 3:3 still breaks, but conflict is counted
	res := phase.Run(ds, cfg)
	if res.Stats.ReadEdgeConflicts == 0 || res.Stats.ReadConflictsAdjudicated == 0 {
		t.Fatalf("read/trio disagreement must be counted: conflicts=%d adjudicated=%d",
			res.Stats.ReadEdgeConflicts, res.Stats.ReadConflictsAdjudicated)
	}
	writeRereadRewrite(t, ds, cfg)
}

// Low-coverage sites retain the low_cov flag through the round trip.
func TestLowCoverageRoundTrip(t *testing.T) {
	ds := scenarioDataset()
	ds.Reads = conflictingReads()
	cfg := config.Default()
	cfg.MinSiteDepth = 10
	res := phase.Run(ds, cfg)
	flagged := 0
	for _, sr := range []string{"child", "father", "mother"} {
		for _, sp := range res.Samples[sr].Sites {
			if sp.LowCov {
				flagged++
			}
		}
	}
	if flagged == 0 {
		t.Fatal("min_site_depth=10 must flag sites as low coverage")
	}
	writeRereadRewrite(t, ds, cfg)
}

func assertRowStatus(t *testing.T, content []byte, sample string, idx int, status, violation string) {
	t.Helper()
	prefix := sample + "\t" + strconv.Itoa(idx) + "\t"
	for _, line := range splitLines(content) {
		if len(line) >= len(prefix) && string(line[:len(prefix)]) == prefix {
			fields := bytes.Split(line, []byte("\t"))
			if string(fields[3]) != status || string(fields[9]) != violation {
				t.Fatalf("site %s/%d: got status=%s violation=%s, want %s/%s; row=%s",
					sample, idx, fields[3], fields[9], status, violation, line)
			}
			return
		}
	}
	t.Fatalf("phase row not found for %s/%d", sample, idx)
}

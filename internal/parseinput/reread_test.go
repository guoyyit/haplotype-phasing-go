package parseinput

import (
	"bytes"
	"encoding/json"
	"strings"
	"testing"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/output"
	"github.com/guoyyit/haplotype-phasing-go/internal/phase"
)

// roundTrip runs the full pipeline (parse -> phase -> write -> re-read ->
// re-phase -> re-write) and returns the original and regenerated phased
// outputs. The two must be byte-identical for every supported input shape.
func roundTrip(t *testing.T, variants, genotypes, reads, trios, cfgJSON string) (orig, regen string) {
	t.Helper()
	v := writeTmp(t, "v.tsv", variants)
	g := writeTmp(t, "g.tsv", genotypes)
	r := ""
	if reads != "" {
		r = writeTmp(t, "r.tsv", reads)
	}
	tr := ""
	if trios != "" {
		tr = writeTmp(t, "t.tsv", trios)
	}
	ds, _, err := LoadDataset(v, g, r, tr)
	if err != nil {
		t.Fatal(err)
	}
	cfg := config.Default()
	if cfgJSON != "" {
		if err := json.Unmarshal([]byte(cfgJSON), &cfg); err != nil {
			t.Fatal(err)
		}
	}
	res := phase.Run(ds, cfg)
	var buf bytes.Buffer
	if err := output.WritePhased(&buf, ds, res, cfg); err != nil {
		t.Fatal(err)
	}
	orig = buf.String()
	phasedPath := writeTmp(t, "phased.tsv", orig)
	ds2, _, cfg2, _, err := LoadPhased(phasedPath)
	if err != nil {
		t.Fatal(err)
	}
	var buf2 bytes.Buffer
	if err := output.WritePhased(&buf2, ds2, phase.Run(ds2, cfg2), cfg2); err != nil {
		t.Fatal(err)
	}
	regen = buf2.String()
	return orig, regen
}

// phaseRows returns the #phase data rows (excluding the header) keyed by
// "sample\tvariant_index".
func phaseRows(t *testing.T, phased string) map[string]string {
	t.Helper()
	rows := map[string]string{}
	inPhase := false
	for _, line := range strings.Split(phased, "\n") {
		line = strings.TrimRight(line, "\r")
		if strings.HasPrefix(line, "#phase") {
			inPhase = true
			continue
		}
		if inPhase && strings.HasPrefix(line, "#") {
			break
		}
		if inPhase && line != "" && !strings.HasPrefix(line, "sample\t") {
			fields := strings.Split(line, "\t")
			rows[fields[0]+"\t"+fields[1]] = line
		}
	}
	return rows
}

// statsSection returns the #stats key/value pairs.
func statsSection(t *testing.T, phased string) map[string]string {
	t.Helper()
	stats := map[string]string{}
	inStats := false
	for _, line := range strings.Split(phased, "\n") {
		line = strings.TrimRight(line, "\r")
		if strings.HasPrefix(line, "#stats") {
			inStats = true
			continue
		}
		if inStats && line != "" {
			fields := strings.Split(line, "\t")
			if len(fields) == 2 {
				stats[fields[0]] = fields[1]
			}
		}
	}
	return stats
}

const variants2 = `#variants
v1	chr1	100	A	G
v2	chr1	200	C	T
`

const trios2 = `#trios
# child	father	mother
child	father	mother
`

// Violated heterozygous and homozygous sites must round-trip byte-identical:
// the violation marker, status and all counts are unchanged across re-read.
func TestRereadViolatedHetAndHomo(t *testing.T) {
	genotypes := `#genotypes
child	v1	A/G
child	v2	C/C
father	v1	G/G
father	v2	T/T
mother	v1	G/G
mother	v2	T/T
`
	orig, regen := roundTrip(t, variants2, genotypes, "", trios2, "")
	if orig != regen {
		t.Fatalf("violated sites do not round-trip byte-identical\n--- orig ---\n%s\n--- regen ---\n%s", orig, regen)
	}
	rows := phaseRows(t, orig)
	// v1: violated het -> UNPHASED, h1=A0, h2=A1, violation=1.
	het := rows["child\t0"]
	if !strings.Contains(het, "\tUNPHASED\t") || !strings.HasSuffix(het, "\t1") {
		t.Fatalf("violated het row wrong: %s", het)
	}
	if !strings.Contains(het, "\tA0\tA1\t") {
		t.Fatalf("violated het must record canonical allele calls: %s", het)
	}
	// v2: violated homo -> UNPHASED, h1=A0, h2=A0, violation=1.
	homo := rows["child\t1"]
	if !strings.Contains(homo, "\tUNPHASED\t") || !strings.HasSuffix(homo, "\t1") {
		t.Fatalf("violated homo row wrong: %s", homo)
	}
	if !strings.Contains(homo, "\tA0\tA0\t") {
		t.Fatalf("violated homo must record canonical allele calls: %s", homo)
	}
	stats := statsSection(t, orig)
	if stats["trio_violations"] != "2" {
		t.Fatalf("expected 2 violations, got %s", stats["trio_violations"])
	}
	if stats["missing_genotypes"] != "0" {
		t.Fatalf("violated sites must not be counted as missing, got %s", stats["missing_genotypes"])
	}
}

// Missing notation on either side (./A and A/.) must parse as partial and
// output PARTIAL with the observed allele preserved.
func TestRereadMissingNotationBothSides(t *testing.T) {
	genotypes := `#genotypes
s1	v1	./A
s1	v2	T/.
`
	orig, regen := roundTrip(t, variants2, genotypes, "", "", "")
	if orig != regen {
		t.Fatalf("partial sites do not round-trip byte-identical\n--- orig ---\n%s\n--- regen ---\n%s", orig, regen)
	}
	rows := phaseRows(t, orig)
	// ./A at v1 (ref=A, alt=G): observed A -> PARTIAL, h1=A0, h2=., a0=A, a1=G.
	left := rows["s1\t0"]
	if !strings.Contains(left, "\tPARTIAL\t") {
		t.Fatalf("./A must be PARTIAL: %s", left)
	}
	if !strings.Contains(left, "\tA\tG\tA0\t.\t") {
		t.Fatalf("./A must preserve ref/alt and observed allele: %s", left)
	}
	// T/. at v2 (ref=C, alt=T): observed T -> PARTIAL, h1=A1, h2=., a0=C, a1=T.
	right := rows["s1\t1"]
	if !strings.Contains(right, "\tPARTIAL\t") {
		t.Fatalf("T/. must be PARTIAL: %s", right)
	}
	if !strings.Contains(right, "\tC\tT\tA1\t.\t") {
		t.Fatalf("T/. must preserve ref/alt and observed allele: %s", right)
	}
}

// Fully missing genotype (./.) must output MISSING and round-trip.
func TestRereadFullyMissing(t *testing.T) {
	genotypes := `#genotypes
s1	v1	./.
`
	orig, regen := roundTrip(t, variants2, genotypes, "", "", "")
	if orig != regen {
		t.Fatalf("fully missing site does not round-trip byte-identical\n--- orig ---\n%s\n--- regen ---\n%s", orig, regen)
	}
	rows := phaseRows(t, orig)
	row := rows["s1\t0"]
	if !strings.Contains(row, "\tMISSING\t") {
		t.Fatalf("./. must be MISSING: %s", row)
	}
	if !strings.Contains(row, "\t.\t.\t") {
		t.Fatalf("MISSING must have no allele calls: %s", row)
	}
	stats := statsSection(t, orig)
	// v1 is explicitly ./. and v2 is absent from the genotype file, so both
	// are fully missing.
	if stats["missing_genotypes"] != "2" {
		t.Fatalf("expected 2 missing genotypes, got %s", stats["missing_genotypes"])
	}
}

// Low-coverage sites must preserve the low_cov flag and round-trip.
func TestRereadLowCoverage(t *testing.T) {
	genotypes := `#genotypes
s1	v1	A/G
s1	v2	C/T
`
	cfgJSON := `{"min_edge_support":2,"min_edge_strength":0.8,"min_block_confidence":0.95,"block_max_error_rate":0.25,"min_site_depth":5,"seed":1,"enable_trio":true,"enable_reads":true}`
	orig, regen := roundTrip(t, variants2, genotypes, "", "", cfgJSON)
	if orig != regen {
		t.Fatalf("low-coverage sites do not round-trip byte-identical\n--- orig ---\n%s\n--- regen ---\n%s", orig, regen)
	}
	rows := phaseRows(t, orig)
	for _, key := range []string{"s1\t0", "s1\t1"} {
		if !strings.Contains(rows[key], "\t1\t") {
			t.Fatalf("site %s must be flagged low coverage: %s", key, rows[key])
		}
	}
}

// Non-default config and seed must round-trip byte-identical.
func TestRereadNonDefaultConfigAndSeed(t *testing.T) {
	genotypes := `#genotypes
s1	v1	A/G
s1	v2	C/T
`
	reads := `#reads
r1	s1	v1	A
r1	s1	v2	T
r2	s1	v1	G
r2	s1	v2	C
`
	cfgJSON := `{"min_edge_support":1,"min_edge_strength":0.6,"min_block_confidence":0.9,"block_max_error_rate":0.4,"min_site_depth":1,"seed":42,"enable_trio":true,"enable_reads":true}`
	orig, regen := roundTrip(t, variants2, genotypes, reads, "", cfgJSON)
	if orig != regen {
		t.Fatalf("non-default config/seed do not round-trip byte-identical\n--- orig ---\n%s\n--- regen ---\n%s", orig, regen)
	}
}

// Disabling reads or trio evidence must still round-trip byte-identical.
func TestRereadReadsOffAndTriosOff(t *testing.T) {
	genotypes := `#genotypes
child	v1	A/G
child	v2	C/T
father	v1	G/G
father	v2	T/T
mother	v1	G/G
mother	v2	C/C
`
	reads := `#reads
r1	child	v1	A
r1	child	v2	T
`
	cfgJSON := `{"min_edge_support":2,"min_edge_strength":0.8,"min_block_confidence":0.95,"block_max_error_rate":0.25,"min_site_depth":3,"seed":1,"enable_trio":false,"enable_reads":false}`
	orig, regen := roundTrip(t, variants2, genotypes, reads, trios2, cfgJSON)
	if orig != regen {
		t.Fatalf("reads-off/trios-off do not round-trip byte-identical\n--- orig ---\n%s\n--- regen ---\n%s", orig, regen)
	}
}

// A genotype allele that is not ref/alt must be rejected (treated as
// missing), never silently rewritten and never emitted as a third allele.
func TestRereadThirdAlleleRejected(t *testing.T) {
	genotypes := `#genotypes
s1	v1	A/C
s1	v2	G/G
`
	orig, regen := roundTrip(t, variants2, genotypes, "", "", "")
	if orig != regen {
		t.Fatalf("third-allele genotypes do not round-trip byte-identical\n--- orig ---\n%s\n--- regen ---\n%s", orig, regen)
	}
	rows := phaseRows(t, orig)
	// v1 A/C at ref=A alt=G: C rejected -> A/. -> PARTIAL, observed A.
	one := rows["s1\t0"]
	if !strings.Contains(one, "\tPARTIAL\t") {
		t.Fatalf("A/C must become PARTIAL after rejecting C: %s", one)
	}
	if !strings.Contains(one, "\tA\tG\tA0\t.\t") {
		t.Fatalf("rejected allele must not enter output: %s", one)
	}
	// v2 G/G at ref=C alt=T: both rejected -> ./. -> MISSING.
	two := rows["s1\t1"]
	if !strings.Contains(two, "\tMISSING\t") {
		t.Fatalf("C/C must become MISSING after rejecting both: %s", two)
	}
	if !strings.Contains(two, "\t.\t.\t") {
		t.Fatalf("rejected alleles must not enter output: %s", two)
	}
}

// Every supported input shape must round-trip byte-identical: this is the
// acceptance criterion for the re-read verification mode.
func TestRereadByteIdenticalAcrossScenarios(t *testing.T) {
	scenarios := []struct {
		name      string
		genotypes string
		reads     string
		trios     string
		cfgJSON   string
	}{
		{
			name: "violated het and homo",
			genotypes: `#genotypes
child	v1	A/G
child	v2	C/C
father	v1	G/G
father	v2	T/T
mother	v1	G/G
mother	v2	T/T
`,
			trios: trios2,
		},
		{
			name: "partial missing both sides",
			genotypes: `#genotypes
s1	v1	./A
s1	v2	T/.
`,
		},
		{
			name: "fully missing",
			genotypes: `#genotypes
s1	v1	./.
s1	v2	./.
`,
		},
		{
			name: "low coverage",
			genotypes: `#genotypes
s1	v1	A/G
s1	v2	C/T
`,
			cfgJSON: `{"min_edge_support":2,"min_edge_strength":0.8,"min_block_confidence":0.95,"block_max_error_rate":0.25,"min_site_depth":5,"seed":1,"enable_trio":true,"enable_reads":true}`,
		},
		{
			name: "reads and trio conflict",
			genotypes: `#genotypes
child	v1	A/G
child	v2	C/T
father	v1	A/A
father	v2	T/T
mother	v1	G/G
mother	v2	C/C
`,
			reads: `#reads
r1	child	v1	A
r1	child	v2	C
`,
			trios: trios2,
		},
		{
			name:      "no reads no trios",
			genotypes: `#genotypes
s1	v1	A/G
s1	v2	C/T
`,
		},
		{
			name: "reads off trios off",
			genotypes: `#genotypes
child	v1	A/G
child	v2	C/T
father	v1	G/G
father	v2	T/T
mother	v1	G/G
mother	v2	C/C
`,
			reads: `#reads
r1	child	v1	A
r1	child	v2	T
`,
			trios:   trios2,
			cfgJSON: `{"min_edge_support":2,"min_edge_strength":0.8,"min_block_confidence":0.95,"block_max_error_rate":0.25,"min_site_depth":3,"seed":7,"enable_trio":false,"enable_reads":false}`,
		},
		{
			name: "non-default seed",
			genotypes: `#genotypes
s1	v1	A/G
s1	v2	C/T
`,
			reads: `#reads
r1	s1	v1	A
r1	s1	v2	T
r2	s1	v1	G
r2	s1	v2	C
`,
			cfgJSON: `{"min_edge_support":2,"min_edge_strength":0.8,"min_block_confidence":0.95,"block_max_error_rate":0.25,"min_site_depth":3,"seed":99,"enable_trio":true,"enable_reads":true}`,
		},
	}
	for _, sc := range scenarios {
		t.Run(sc.name, func(t *testing.T) {
			orig, regen := roundTrip(t, variants2, sc.genotypes, sc.reads, sc.trios, sc.cfgJSON)
			if orig != regen {
				t.Fatalf("round-trip produced a different result\n--- orig ---\n%s\n--- regen ---\n%s", orig, regen)
			}
		})
	}
}

// Legacy 1.0 files with violated heterozygous sites (h1/h2 = ".") cannot be
// re-read correctly; the loader must give a clear, locatable compatibility
// error instead of an unlocatable verification failure.
func TestRereadLegacyViolatedHetGivesClearError(t *testing.T) {
	legacy := `#haplophased	1.0
#variants
variant_id	chrom	pos	ref	alt
v1	chr1	100	A	G
#trios
child	father	mother
child	father	mother
#reads
read_id	sample	variant_id	allele
#config
{"min_edge_support":2,"min_edge_strength":0.8,"min_block_confidence":0.95,"block_max_error_rate":0.25,"min_site_depth":3,"seed":1,"enable_trio":true,"enable_reads":true}
#phase
sample	variant_index	variant_id	status	low_cov	a0	a1	h1	h2	violation
child	0	v1	UNPHASED	1	A	G	.	.	1
#blocks
sample	block_id	start_idx	end_idx	site_count	hets_count	support	conflict	confidence	error_rate	status
#stats
samples	3
variants	1
hets_total	1
phased_hets	0
unphased_hets	1
blocks_total	0
read_observations	0
read_pairs	0
read_edge_supports	0
read_edge_conflicts	0
read_conflicts_adjudicated	0
trio_count	1
trio_constraint_sites	1
trio_violations	1
breaks_by_threshold	0
low_coverage_sites	3
missing_genotypes	0
`
	path := writeTmp(t, "legacy.tsv", legacy)
	_, _, _, _, err := LoadPhased(path)
	if err == nil {
		t.Fatal("legacy violated het file must be rejected with a compatibility error")
	}
	msg := err.Error()
	if !strings.Contains(msg, "legacy phased file") || !strings.Contains(msg, "v1") {
		t.Fatalf("error must be locatable and mention the site id: %v", err)
	}
}

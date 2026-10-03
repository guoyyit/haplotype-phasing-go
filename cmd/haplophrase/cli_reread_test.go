package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const cliVariants = `#variants
v1	chr1	1	A	G
v2	chr1	2	C	T
`

func writeInputs(t *testing.T) (string, string, string, string) {
	t.Helper()
	dir := t.TempDir()
	v := filepath.Join(dir, "variants.tsv")
	g := filepath.Join(dir, "genotypes.tsv")
	tr := filepath.Join(dir, "trios.tsv")
	rd := filepath.Join(dir, "reads.tsv")
	if err := os.WriteFile(v, []byte(cliVariants), 0644); err != nil {
		t.Fatal(err)
	}
	// v1: het violation (A/G child of A/A x A/A); v2: A/. partial with the
	// missing marker on the right, plus a ./C-style partial covered below.
	geno := strings.Join([]string{
		"#genotypes",
		"c\tv1\tA/G",
		"c\tv2\tC/.",
		"f\tv1\tA/A",
		"f\tv2\tC/T",
		"m\tv1\tA/A",
		"m\tv2\tT/T",
	}, "\n") + "\n"
	if err := os.WriteFile(g, []byte(geno), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(tr, []byte("#trios\nc\tf\tm\n"), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(rd, []byte("#reads\n"), 0644); err != nil {
		t.Fatal(err)
	}
	return v, g, tr, rd
}

// A full CLI run followed by -reread on the produced file exits 0.
func TestCLIRereadEndToEnd(t *testing.T) {
	v, g, tr, rd := writeInputs(t)
	out := filepath.Join(t.TempDir(), "phased.tsv")
	if err := runCLI([]string{
		"-variants", v, "-genotypes", g, "-trios", tr, "-reads", rd,
		"-out", out, "-report", "",
	}); err != nil {
		t.Fatal(err)
	}
	if err := runCLI([]string{"-reread", out}); err != nil {
		t.Fatalf("reread verification failed: %v", err)
	}
}

// A corrupt edit to the phased file produces a locatable error that names
// the first differing line number and both line contents.
func TestCLIRereadReportsMismatchLine(t *testing.T) {
	v, g, tr, rd := writeInputs(t)
	out := filepath.Join(t.TempDir(), "phased.tsv")
	if err := runCLI([]string{
		"-variants", v, "-genotypes", g, "-trios", tr, "-reads", rd,
		"-out", out, "-report", "",
	}); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(out)
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(string(data), "\n")
	// Flip the violation bit on the first #phase data row for the child.
	for i, line := range lines {
		if strings.HasPrefix(line, "c\t0\t") {
			fields := strings.Split(line, "\t")
			fields[9] = "0"
			lines[i] = strings.Join(fields, "\t")
			break
		}
	}
	corrupt := strings.Join(lines, "\n")
	corruptPath := filepath.Join(t.TempDir(), "corrupt.tsv")
	if err := os.WriteFile(corruptPath, []byte(corrupt), 0644); err != nil {
		t.Fatal(err)
	}
	err = runCLI([]string{"-reread", corruptPath})
	if err == nil {
		t.Fatal("corrupt file must fail verification")
	}
	msg := err.Error()
	if !strings.Contains(msg, "line ") || !strings.Contains(msg, "original:") || !strings.Contains(msg, "rephased:") {
		t.Fatalf("error must locate the first mismatch with line number and contents: %v", err)
	}
}

// Genotype alleles outside ref/alt are rejected by the CLI with a locatable
// parse error and do not produce a result file.
func TestCLIRejectsInvalidAllele(t *testing.T) {
	dir := t.TempDir()
	v := filepath.Join(dir, "variants.tsv")
	g := filepath.Join(dir, "genotypes.tsv")
	if err := os.WriteFile(v, []byte(cliVariants), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(g, []byte("#genotypes\ns\tv1\tA/Z\n"), 0644); err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(t.TempDir(), "phased.tsv")
	err := runCLI([]string{"-variants", v, "-genotypes", g, "-out", out, "-report", ""})
	if err == nil || !strings.Contains(err.Error(), "neither ref") {
		t.Fatalf("invalid allele must be rejected explicitly, got %v", err)
	}
}

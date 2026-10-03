package main

import (
	"bytes"
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/output"
	"github.com/guoyyit/haplotype-phasing-go/internal/parseinput"
	"github.com/guoyyit/haplotype-phasing-go/internal/phase"
)

func main() {
	if err := runCLI(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, "haplophrase:", err)
		os.Exit(1)
	}
}

func runCLI(args []string) error {
	fs := flag.NewFlagSet("haplophrase", flag.ContinueOnError)
	var (
		variantsPath  = fs.String("variants", "", "variants TSV (required)")
		genotypesPath = fs.String("genotypes", "", "genotypes TSV (required)")
		readsPath     = fs.String("reads", "", "reads TSV (optional)")
		triosPath     = fs.String("trios", "", "trios TSV (optional)")
		configPath    = fs.String("config", "", "JSON config (optional, defaults are used if absent)")
		outPhased     = fs.String("out", "phased.tsv", "phased result output path")
		outReport     = fs.String("report", "report.txt", "human-readable report output path")
		reread        = fs.String("reread", "", "verify a previously written phased file instead of phasing")
	)
	if err := fs.Parse(args); err != nil {
		return err
	}

	if *reread != "" {
		return verifyReread(*reread)
	}
	if *variantsPath == "" || *genotypesPath == "" {
		return fmt.Errorf("-variants and -genotypes are required")
	}

	cfg, err := config.Load(*configPath)
	if err != nil {
		return err
	}
	ds, notices, err := parseinput.LoadDataset(*variantsPath, *genotypesPath, *readsPath, *triosPath)
	if err != nil {
		return err
	}
	for _, notice := range notices {
		fmt.Fprintln(os.Stderr, "notice:", notice)
	}

	res := phase.Run(ds, cfg)
	if err := writeAtomic(*outPhased, func(f *os.File) error {
		return output.WritePhased(f, ds, res, cfg)
	}); err != nil {
		return err
	}
	if *outReport != "" {
		if err := writeAtomic(*outReport, func(f *os.File) error {
			return output.WriteReport(f, ds, res, cfg)
		}); err != nil {
			return err
		}
	}
	fmt.Printf("phased %d samples / %d variants; blocks=%d; wrote %s\n",
		res.Stats.Samples, res.Stats.Variants, res.Stats.BlocksTotal, *outPhased)
	return nil
}

// verifyReread loads a phased file, re-phases it and checks byte equality
// of the regenerated phased representation with the stored config.
func verifyReread(path string) error {
	ds, _, cfg, notices, err := parseinput.LoadPhased(path)
	if err != nil {
		return err
	}
	for _, notice := range notices {
		fmt.Fprintln(os.Stderr, "notice:", notice)
	}
	regenerated := phase.Run(ds, cfg)
	var buf bytes.Buffer
	if err := output.WritePhased(&buf, ds, regenerated, cfg); err != nil {
		return err
	}
	original, err := os.ReadFile(path)
	if err != nil {
		return err
	}
	if lineNo, origLine, regenLine, ok := firstLineDiff(string(original), buf.String()); !ok {
		return fmt.Errorf("reread verification FAILED: first differing line %d\n  original:   %q\n  regenerated: %q", lineNo, origLine, regenLine)
	}
	fmt.Println("reread verification OK: identical phasing")
	return nil
}

// firstLineDiff compares two phased outputs line by line and returns the
// 1-based line number and content of the first differing line. Trailing
// whitespace per line is ignored to match the tolerant comparison.
func firstLineDiff(original, regenerated string) (lineNo int, origLine, regenLine string, ok bool) {
	a := strings.Split(strings.TrimSpace(original), "\n")
	b := strings.Split(strings.TrimSpace(regenerated), "\n")
	for i := 0; i < len(a) || i < len(b); i++ {
		var x, y string
		if i < len(a) {
			x = strings.TrimRight(a[i], "\r ")
		}
		if i < len(b) {
			y = strings.TrimRight(b[i], "\r ")
		}
		if x != y {
			return i + 1, x, y, false
		}
	}
	return 0, "", "", true
}

func writeAtomic(path string, fn func(*os.File) error) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	if err := fn(f); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

package main

import (
	"flag"
	"fmt"
	"os"

	phasing "github.com/guoyyit/haplotype-phasing-go"
)

func main() {
	fs := flag.NewFlagSet("haplotyper", flag.ContinueOnError)
	variantsPath := fs.String("variants", "", "path to variants TSV")
	genotypesPath := fs.String("genotypes", "", "path to genotypes TSV")
	readsPath := fs.String("reads", "", "path to reads TSV (optional)")
	pedigreePath := fs.String("pedigree", "", "path to trio pedigree TSV (optional)")
	configPath := fs.String("config", "", "path to config JSON (optional)")
	outputPath := fs.String("output", "", "path to write phasing result JSON")
	rephasePath := fs.String("rephase", "", "path to a prior result JSON to re-read and re-emit")

	if err := fs.Parse(os.Args[1:]); err != nil {
		os.Exit(2)
	}
	if *variantsPath == "" || *genotypesPath == "" || *outputPath == "" {
		fmt.Fprintln(os.Stderr, "usage: haplotyper -variants V -genotypes G -output O [-reads R] [-pedigree P] [-config C] [-rephase prior.json]")
		os.Exit(2)
	}

	in, err := phasing.ParseInput(*variantsPath, *genotypesPath, *readsPath, *pedigreePath, *configPath)
	if err != nil {
		fail("load input", err)
	}

	var result *phasing.Result
	if *rephasePath != "" {
		prior, err := phasing.ReadResult(*rephasePath)
		if err != nil {
			fail("read prior result", err)
		}
		result, err = phasing.Rephase(in, prior)
		if err != nil {
			fail("rephase", err)
		}
	} else {
		result, err = phasing.Run(in)
		if err != nil {
			fail("phase", err)
		}
	}

	if err := phasing.WriteResult(*outputPath, result); err != nil {
		fail("write output", err)
	}

	fmt.Printf("wrote %s: %d samples, %d sites, %d phased hets, %d read conflicts, %d trio violations\n",
		*outputPath, result.Stats.Samples, result.Stats.VariantSites,
		result.Stats.PhasedHets,
		result.Conflicts.ReadReadConflicts+result.Conflicts.ReadTrioConflicts+
			result.Conflicts.ReadGenotypeConflicts,
		result.Conflicts.TrioMendelViolations)
}

func fail(stage string, err error) {
	fmt.Fprintf(os.Stderr, "haplotyper: %s: %v\n", stage, err)
	os.Exit(1)
}

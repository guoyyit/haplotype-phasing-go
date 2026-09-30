package phasing

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// ParseInput reads variants, genotypes and optional reads, pedigree and config
// files from the local filesystem. Empty optional paths mean "not provided".
func ParseInput(variantsPath, genotypesPath, readsPath, pedigreePath, configPath string) (*Input, error) {
	in := &Input{}

	cfg := DefaultConfig()
	if configPath != "" {
		raw, err := os.ReadFile(configPath)
		if err != nil {
			return nil, fmt.Errorf("read config: %w", err)
		}
		if err := json.Unmarshal(raw, &cfg); err != nil {
			return nil, fmt.Errorf("parse config: %w", err)
		}
	}
	resolved, err := cfg.Resolve()
	if err != nil {
		return nil, err
	}
	in.Config = resolved.Config

	variants, err := readVariants(variantsPath)
	if err != nil {
		return nil, err
	}
	// Deterministic canonical site order: chromosome, position, id.
	sort.SliceStable(variants, func(i, j int) bool {
		if variants[i].Chrom != variants[j].Chrom {
			return variants[i].Chrom < variants[j].Chrom
		}
		if variants[i].Pos != variants[j].Pos {
			return variants[i].Pos < variants[j].Pos
		}
		return variants[i].ID < variants[j].ID
	})
	for i := range variants {
		variants[i].Index = i
	}
	in.Variants = variants

	samples, err := readGenotypes(genotypesPath, variants)
	if err != nil {
		return nil, err
	}
	in.Samples = samples

	if readsPath != "" {
		reads, err := readReads(readsPath, variants, samples)
		if err != nil {
			return nil, err
		}
		in.Reads = reads
	}

	if pedigreePath != "" {
		trio, err := readTrio(pedigreePath, samples)
		if err != nil {
			return nil, err
		}
		in.Trio = trio
	}

	return in, nil
}

// tsvScanner streams a tab separated file one content row at a time. Blank
// lines and '#' comments are skipped. Streaming keeps peak memory independent
// of file size at the documented scale boundary.
type tsvScanner struct {
	f       *os.File
	scanner *bufio.Scanner
	lineNo  int
	cols    []string
}

func openTSV(path string) (*tsvScanner, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	s := newLineScanner(f)
	return &tsvScanner{f: f, scanner: s}, nil
}

// next advances to the next non-comment, non-blank row.
func (t *tsvScanner) next() bool {
	for t.scanner.Scan() {
		t.lineNo++
		line := strings.TrimRight(t.scanner.Text(), "\r")
		if strings.TrimSpace(line) == "" || strings.HasPrefix(strings.TrimSpace(line), "#") {
			continue
		}
		t.cols = strings.Split(line, "\t")
		return true
	}
	return false
}

func (t *tsvScanner) close() error {
	err := t.scanner.Err()
	if cerr := t.f.Close(); err == nil {
		err = cerr
	}
	return err
}

func readVariants(path string) ([]Variant, error) {
	ts, err := openTSV(path)
	if err != nil {
		return nil, fmt.Errorf("read variants: %w", err)
	}
	defer ts.close()

	if !ts.next() {
		return nil, fmt.Errorf("variants file %s is empty", path)
	}
	header := ts.cols
	if len(header) < 5 || header[0] != "variant_id" || header[1] != "chrom" ||
		header[2] != "pos" || header[3] != "ref" || header[4] != "alt" {
		return nil, fmt.Errorf("variants file %s: expected header variant_id chrom pos ref alt", path)
	}

	seen := make(map[string]bool)
	variants := []Variant{}
	for ts.next() {
		row := ts.cols
		if len(row) < 5 {
			return nil, fmt.Errorf("variants file %s line %d: expected 5 columns", path, ts.lineNo)
		}
		id := row[0]
		if id == "" {
			return nil, fmt.Errorf("variants file %s line %d: empty variant_id", path, ts.lineNo)
		}
		if seen[id] {
			return nil, fmt.Errorf("variants file %s: duplicate variant_id %q", path, id)
		}
		seen[id] = true
		pos, err := strconv.ParseInt(row[2], 10, 64)
		if err != nil || pos < 0 {
			return nil, fmt.Errorf("variants file %s line %d: invalid pos %q", path, ts.lineNo, row[2])
		}
		if row[3] == "" || row[4] == "" {
			return nil, fmt.Errorf("variants file %s line %d: ref/alt must be non-empty", path, ts.lineNo)
		}
		variants = append(variants, Variant{ID: id, Chrom: row[1], Pos: pos, Ref: row[3], Alt: row[4]})
	}
	if err := ts.close(); err != nil {
		return nil, fmt.Errorf("read variants: %w", err)
	}
	return variants, nil
}

func readGenotypes(path string, variants []Variant) ([]Sample, error) {
	ts, err := openTSV(path)
	if err != nil {
		return nil, fmt.Errorf("read genotypes: %w", err)
	}
	defer ts.close()

	if !ts.next() {
		return nil, fmt.Errorf("genotypes file %s is empty", path)
	}
	header := ts.cols
	if len(header) < 1 || header[0] != "sample_id" {
		return nil, fmt.Errorf("genotypes file %s: first header column must be sample_id", path)
	}
	idxByID := make(map[string]int, len(variants))
	for _, v := range variants {
		idxByID[v.ID] = v.Index
	}
	colIdx := make([]int, len(header)-1)
	for c := 1; c < len(header); c++ {
		idx, ok := idxByID[header[c]]
		if !ok {
			return nil, fmt.Errorf("genotypes file %s: unknown variant_id %q in header", path, header[c])
		}
		colIdx[c-1] = idx
	}

	samples := []Sample{}
	seenSamples := make(map[string]bool)
	for ts.next() {
		row := ts.cols
		if len(row)-1 != len(colIdx) {
			return nil, fmt.Errorf("genotypes file %s line %d: expected %d genotype columns, got %d",
				path, ts.lineNo, len(colIdx), len(row)-1)
		}
		id := row[0]
		if id == "" {
			return nil, fmt.Errorf("genotypes file %s line %d: empty sample_id", path, ts.lineNo)
		}
		if seenSamples[id] {
			return nil, fmt.Errorf("genotypes file %s: duplicate sample_id %q", path, id)
		}
		seenSamples[id] = true
		gt := make([]GenotypeCall, len(variants))
		for c, vi := range colIdx {
			call, err := parseGenotypeToken(row[c+1])
			if err != nil {
				return nil, fmt.Errorf("genotypes file %s line %d sample %q variant %q: %w",
					path, ts.lineNo, id, variants[vi].ID, err)
			}
			gt[vi] = call
		}
		samples = append(samples, Sample{ID: id, GT: gt})
	}
	if err := ts.close(); err != nil {
		return nil, fmt.Errorf("read genotypes: %w", err)
	}
	return samples, nil
}

// parseGenotypeToken parses tokens such as "0/1", "./1" and "1|1".
func parseGenotypeToken(tok string) (GenotypeCall, error) {
	if tok == "." {
		return GenotypeCall{Present: true, A1: -1, A2: -1}, nil
	}
	if strings.Count(tok, "/")+strings.Count(tok, "|") != 1 {
		return GenotypeCall{}, fmt.Errorf("invalid genotype %q (expected A/B)", tok)
	}
	sep := "/"
	if strings.Contains(tok, "|") {
		sep = "|"
	}
	parts := strings.Split(tok, sep)
	a1, err1 := parseAlleleToken(parts[0])
	a2, err2 := parseAlleleToken(parts[1])
	if err1 != nil || err2 != nil {
		return GenotypeCall{}, fmt.Errorf("invalid genotype %q", tok)
	}
	return GenotypeCall{Present: true, A1: a1, A2: a2}, nil
}

func parseAlleleToken(s string) (int8, error) {
	switch s {
	case "0":
		return 0, nil
	case "1":
		return 1, nil
	case ".":
		return -1, nil
	default:
		return 0, fmt.Errorf("allele must be 0, 1 or ., got %q", s)
	}
}

func readReads(path string, variants []Variant, samples []Sample) ([]Read, error) {
	ts, err := openTSV(path)
	if err != nil {
		return nil, fmt.Errorf("read reads: %w", err)
	}
	defer ts.close()

	if !ts.next() {
		return nil, fmt.Errorf("reads file %s is empty", path)
	}
	header := ts.cols
	if len(header) < 4 || header[0] != "read_id" || header[1] != "sample_id" ||
		header[2] != "variant_id" || header[3] != "allele" {
		return nil, fmt.Errorf("reads file %s: expected header read_id sample_id variant_id allele [weight]", path)
	}
	sampleIdx := make(map[string]int, len(samples))
	for i, s := range samples {
		sampleIdx[s.ID] = i
	}
	variantIdx := make(map[string]int, len(variants))
	for _, v := range variants {
		variantIdx[v.ID] = v.Index
	}

	// Stream into compact records first; records carry a sample ordinal rather
	// than holding sample strings, and read identities are interned to ids.
	type rec struct {
		readID   int
		sampleID int
		variant  int
		allele   int8
		weight   float64
	}
	readNames := []string{}
	readIntern := map[string]int{}
	readSample := map[int]int{}
	var recs []rec

	for ts.next() {
		row := ts.cols
		if len(row) < 4 {
			return nil, fmt.Errorf("reads file %s line %d: expected >= 4 columns", path, ts.lineNo)
		}
		readIDStr, sampleIDStr, variantIDStr := row[0], row[1], row[2]
		si, ok := sampleIdx[sampleIDStr]
		if !ok {
			return nil, fmt.Errorf("reads file %s line %d: unknown sample_id %q", path, ts.lineNo, sampleIDStr)
		}
		vi, ok := variantIdx[variantIDStr]
		if !ok {
			return nil, fmt.Errorf("reads file %s line %d: unknown variant_id %q", path, ts.lineNo, variantIDStr)
		}
		al, err := parseAlleleToken(row[3])
		if err != nil || al < 0 {
			return nil, fmt.Errorf("reads file %s line %d: allele must be 0 or 1, got %q", path, ts.lineNo, row[3])
		}
		weight := 1.0
		if len(row) >= 5 && row[4] != "" {
			weight, err = strconv.ParseFloat(row[4], 64)
			if err != nil || weight <= 0 {
				return nil, fmt.Errorf("reads file %s line %d: invalid weight %q", path, ts.lineNo, row[4])
			}
		}
		rid, seen := readIntern[readIDStr]
		if !seen {
			rid = len(readNames)
			readIntern[readIDStr] = rid
			readNames = append(readNames, readIDStr)
			readSample[rid] = si
		} else if readSample[rid] != si {
			return nil, fmt.Errorf("reads file %s line %d: read %q assigned to multiple samples",
				path, ts.lineNo, readIDStr)
		}
		recs = append(recs, rec{readID: rid, sampleID: si, variant: vi, allele: al, weight: weight})
	}
	if err := ts.close(); err != nil {
		return nil, fmt.Errorf("read reads: %w", err)
	}

	// Canonical, stable order: (sample ordinal, read id, variant). Duplicate
	// (read, variant) rows are validated and merged here; after sorting the
	// records are streamed into Read values without retaining intermediate maps.
	sort.Slice(recs, func(i, j int) bool {
		if recs[i].sampleID != recs[j].sampleID {
			return recs[i].sampleID < recs[j].sampleID
		}
		if recs[i].readID != recs[j].readID {
			return recs[i].readID < recs[j].readID
		}
		return recs[i].variant < recs[j].variant
	})

	reads := []Read{}
	var cur *Read
	flush := func() {
		if cur != nil {
			reads = append(reads, *cur)
			cur = nil
		}
	}
	for i := range recs {
		r := &recs[i]
		if cur == nil || cur.Sample != samples[r.sampleID].ID || cur.ID != readNames[r.readID] {
			flush()
			cur = &Read{ID: readNames[r.readID], Sample: samples[r.sampleID].ID}
		}
		if n := len(cur.Obs); n > 0 && cur.Obs[n-1].Variant == r.variant {
			prev := &cur.Obs[n-1]
			if prev.Allele != r.allele {
				return nil, fmt.Errorf("reads file: read %q reports conflicting alleles for variant %q",
					cur.ID, variants[r.variant].ID)
			}
			prev.Weight += r.weight
			continue
		}
		cur.Obs = append(cur.Obs, SiteObs{Variant: r.variant, Allele: r.allele, Weight: r.weight})
	}
	flush()
	return reads, nil
}

func readTrio(path string, samples []Sample) (*Trio, error) {
	ts, err := openTSV(path)
	if err != nil {
		return nil, fmt.Errorf("read pedigree: %w", err)
	}
	defer ts.close()

	if !ts.next() {
		return nil, fmt.Errorf("pedigree file %s must contain a header and one trio row", path)
	}
	header := ts.cols
	if len(header) < 3 || header[0] != "child" || header[1] != "father" || header[2] != "mother" {
		return nil, fmt.Errorf("pedigree file %s: expected header child father mother", path)
	}
	if !ts.next() {
		return nil, fmt.Errorf("pedigree file %s must contain one trio row", path)
	}
	row := ts.cols
	if ts.next() {
		return nil, fmt.Errorf("pedigree file %s: exactly one trio row is supported", path)
	}
	known := make(map[string]bool, len(samples))
	for _, s := range samples {
		known[s.ID] = true
	}
	trio := &Trio{Child: row[0], Father: row[1], Mother: row[2]}
	roles := []struct{ role, id string }{
		{"child", trio.Child},
		{"father", trio.Father},
		{"mother", trio.Mother},
	}
	for _, member := range roles {
		role, id := member.role, member.id
		if id == "" {
			return nil, fmt.Errorf("pedigree file %s: empty %s sample id", path, role)
		}
		if !known[id] {
			return nil, fmt.Errorf("pedigree file %s: %s sample %q not present in genotypes", path, role, id)
		}
	}
	return trio, nil
}

// newLineScanner returns a scanner with a large buffer so long TSV rows
// (many samples or variants) do not trigger bufio.Scanner token errors.
func newLineScanner(r io.Reader) *bufio.Scanner {
	s := bufio.NewScanner(r)
	s.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	return s
}

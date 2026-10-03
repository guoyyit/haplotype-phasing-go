package parseinput

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

// LoadDataset parses the four input files. Empty paths mean "not provided".
// Notices are deterministic, non-fatal problems (bad allele, unknown site).
func LoadDataset(variantsPath, genotypesPath, readsPath, triosPath string) (ds *model.Dataset, notices []string, err error) {
	ds = &model.Dataset{GT: map[string][]model.Genotype{}}
	if err = parseVariants(variantsPath, ds); err != nil {
		return nil, nil, err
	}
	if err = parseGenotypes(genotypesPath, ds, &notices); err != nil {
		return nil, nil, err
	}
	if readsPath != "" {
		if err = parseReads(readsPath, ds, &notices); err != nil {
			return nil, nil, err
		}
	}
	if triosPath != "" {
		if err = parseTrios(triosPath, ds); err != nil {
			return nil, nil, err
		}
	}
	return ds, notices, nil
}

// LoadPhased reads a previously written phased result file.
func LoadPhased(path string) (*model.Dataset, *model.Result, config.Config, []string, error) {
	return parsePhased(path)
}

// openLines opens path and returns a line scanner plus the file handle.
func openLines(path string) (*os.File, *bufio.Scanner, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	return file, bufio.NewScanner(file), nil
}

// isHeaderRow reports whether a TSV row is a textual column header by
// matching its first cell against any of the known header tokens.
func isHeaderRow(fields []string, tokens ...string) bool {
	if len(fields) == 0 {
		return false
	}
	for _, token := range tokens {
		if fields[0] == token {
			return true
		}
	}
	return false
}

// parseVariants parses the variants TSV:
// #variants header then variant_id chrom pos ref alt.
func parseVariants(path string, ds *model.Dataset) error {
	file, scanner, err := openLines(path)
	if err != nil {
		return err
	}
	defer file.Close()
	seen := map[string]bool{}
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if isHeaderRow(fields, "variant_id") {
			continue
		}
		if len(fields) < 5 {
			return fmt.Errorf("%s:%d: expected 5 columns, got %d", path, lineNo, len(fields))
		}
		id := fields[0]
		if seen[id] {
			return fmt.Errorf("%s:%d: duplicate variant id %q", path, lineNo, id)
		}
		seen[id] = true
		pos, err := strconv.Atoi(fields[2])
		if err != nil || pos < 0 {
			return fmt.Errorf("%s:%d: invalid position %q", path, lineNo, fields[2])
		}
		ds.Variants = append(ds.Variants, model.Variant{
			ID:    id,
			Chrom: fields[1],
			Pos:   pos,
			Ref:   fields[3],
			Alt:   fields[4],
		})
	}
	return scanner.Err()
}

// parseGenotypes parses: sample_id variant_index gt where gt uses / or |.
func parseGenotypes(path string, ds *model.Dataset, notices *[]string) error {
	file, scanner, err := openLines(path)
	if err != nil {
		return err
	}
	defer file.Close()
	indexOf := map[string]int{}
	for idx, variant := range ds.Variants {
		indexOf[variant.ID] = idx
	}
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if isHeaderRow(fields, "sample_id") {
			continue
		}
		if len(fields) < 3 {
			return fmt.Errorf("%s:%d: expected 3 columns, got %d", path, lineNo, len(fields))
		}
		sampleID, variantID := fields[0], fields[1]
		idx, ok := indexOf[variantID]
		if !ok {
			return fmt.Errorf("%s:%d: unknown variant id %q", path, lineNo, variantID)
		}
		if _, ok := ds.GT[sampleID]; !ok {
			ds.Samples = append(ds.Samples, sampleID)
			ds.GT[sampleID] = make([]model.Genotype, len(ds.Variants))
		}
		gt := parseGenotype(fields[2])
		if !gt.Present {
			*notices = append(*notices, fmt.Sprintf("missing genotype: sample=%s variant=%s line=%d", sampleID, variantID, lineNo))
		}
		if gt.Present {
			gt = validateGenotypeAlleles(gt, ds.Variants[idx], sampleID, variantID, lineNo, notices)
		}
		ds.GT[sampleID][idx] = gt
	}
	return scanner.Err()
}

// validateGenotypeAlleles rejects any allele that is not the variant's
// reference or alternate allele: such an allele is treated as missing
// (rejected, never silently rewritten and never emitted as a third allele)
// and the rejection is reported as a notice. The returned genotype has the
// same ploidy semantics as parseGenotype.
func validateGenotypeAlleles(gt model.Genotype, variant model.Variant, sampleID, variantID string, lineNo int, notices *[]string) model.Genotype {
	if gt.Ploidy == 2 {
		if gt.A1 != variant.Ref && gt.A1 != variant.Alt {
			*notices = append(*notices, fmt.Sprintf("genotype allele rejected (not ref/alt): sample=%s variant=%s allele=%s line=%d", sampleID, variantID, gt.A1, lineNo))
			gt.A1 = ""
		}
		if gt.A2 != variant.Ref && gt.A2 != variant.Alt {
			*notices = append(*notices, fmt.Sprintf("genotype allele rejected (not ref/alt): sample=%s variant=%s allele=%s line=%d", sampleID, variantID, gt.A2, lineNo))
			gt.A2 = ""
		}
		switch {
		case gt.A1 == "" && gt.A2 == "":
			return model.Genotype{Present: true, Ploidy: 0}
		case gt.A1 == "":
			return model.Genotype{Present: true, Ploidy: 1, A1: gt.A2}
		case gt.A2 == "":
			return model.Genotype{Present: true, Ploidy: 1, A1: gt.A1}
		}
		return gt
	}
	if gt.Ploidy == 1 {
		if gt.A1 != variant.Ref && gt.A1 != variant.Alt {
			*notices = append(*notices, fmt.Sprintf("genotype allele rejected (not ref/alt): sample=%s variant=%s allele=%s line=%d", sampleID, variantID, gt.A1, lineNo))
			return model.Genotype{Present: true, Ploidy: 0}
		}
	}
	return gt
}

// parseGenotype accepts ./. , ./A , A/. , A/A , A/B with either / or |
// separator. The missing notation "." (or an empty side) is handled on
// either side: a single observed allele yields ploidy 1, two missing sides
// yield ploidy 0.
func parseGenotype(token string) model.Genotype {
	token = strings.ReplaceAll(token, "|", "/")
	parts := strings.Split(token, "/")
	gt := model.Genotype{}
	if len(parts) != 2 {
		return gt
	}
	a1, a2 := strings.TrimSpace(parts[0]), strings.TrimSpace(parts[1])
	m1 := a1 == "." || a1 == ""
	m2 := a2 == "." || a2 == ""
	switch {
	case m1 && m2:
		return model.Genotype{Present: true, Ploidy: 0}
	case m1:
		return model.Genotype{Present: true, Ploidy: 1, A1: a2}
	case m2:
		return model.Genotype{Present: true, Ploidy: 1, A1: a1}
	}
	return model.Genotype{Present: true, Ploidy: 2, A1: a1, A2: a2}
}

// readKey groups observations of one (read, sample).
type readKey struct {
	sample string
	readID string
}

// parseReads parses: read_id sample_id variant_id allele.
func parseReads(path string, ds *model.Dataset, notices *[]string) error {
	file, scanner, err := openLines(path)
	if err != nil {
		return err
	}
	defer file.Close()
	indexOf := map[string]int{}
	for idx, variant := range ds.Variants {
		indexOf[variant.ID] = idx
	}
	byRead := map[readKey][]model.ReadObservation{}
	var order []readKey
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if isHeaderRow(fields, "read_id") {
			continue
		}
		if len(fields) < 4 {
			return fmt.Errorf("%s:%d: expected 4 columns, got %d", path, lineNo, len(fields))
		}
		key := readKey{sample: fields[1], readID: fields[0]}
		idx, ok := indexOf[fields[2]]
		if !ok {
			*notices = append(*notices, fmt.Sprintf("read observation ignored (unknown variant): read=%s variant=%s line=%d", key.readID, fields[2], lineNo))
			continue
		}
		allele := fields[3]
		if allele == "." || allele == "" {
			*notices = append(*notices, fmt.Sprintf("read observation ignored (missing allele): read=%s variant=%s line=%d", key.readID, fields[2], lineNo))
			continue
		}
		variant := ds.Variants[idx]
		if allele != variant.Ref && allele != variant.Alt {
			*notices = append(*notices, fmt.Sprintf("read observation ignored (allele not ref/alt): read=%s variant=%s allele=%s line=%d", key.readID, fields[2], allele, lineNo))
			continue
		}
		if _, ok := byRead[key]; !ok {
			order = append(order, key)
		}
		byRead[key] = append(byRead[key], model.ReadObservation{VariantIdx: idx, Allele: allele})
	}
	if err := scanner.Err(); err != nil {
		return err
	}
	for _, key := range order {
		obs := byRead[key]
		sort.SliceStable(obs, func(i, j int) bool { return obs[i].VariantIdx < obs[j].VariantIdx })
		ds.Reads = append(ds.Reads, model.Read{
			ID:           key.readID,
			Sample:       key.sample,
			Observations: dedupeObservations(obs, key.readID, notices),
		})
	}
	return nil
}

// dedupeObservations keeps a deterministic allele when a read reports one
// site twice: the allele matching the genotype-independent lexical order
// wins, and the collision is reported as a notice and counted downstream.
func dedupeObservations(obs []model.ReadObservation, readID string, notices *[]string) []model.ReadObservation {
	out := obs[:0]
	for i := 0; i < len(obs); {
		j := i + 1
		for j < len(obs) && obs[j].VariantIdx == obs[i].VariantIdx {
			j++
		}
		chosen := obs[i]
		if j-i > 1 {
			alleles := []string{}
			for k := i; k < j; k++ {
				alleles = append(alleles, obs[k].Allele)
			}
			sort.Strings(alleles)
			chosen.Allele = alleles[0]
			*notices = append(*notices, fmt.Sprintf("read internal conflict resolved deterministically: read=%s variant_index=%d alleles=%s", readID, chosen.VariantIdx, strings.Join(alleles, ",")))
		}
		out = append(out, chosen)
		i = j
	}
	return out
}

// parseTrios parses: child_id father_id mother_id.
func parseTrios(path string, ds *model.Dataset) error {
	file, scanner, err := openLines(path)
	if err != nil {
		return err
	}
	defer file.Close()
	seen := map[string]bool{}
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if isHeaderRow(fields, "child_id") {
			continue
		}
		if len(fields) < 3 {
			return fmt.Errorf("%s:%d: expected 3 columns, got %d", path, lineNo, len(fields))
		}
		if seen[fields[0]] {
			return fmt.Errorf("%s:%d: duplicate trio child %q", path, lineNo, fields[0])
		}
		seen[fields[0]] = true
		ds.Trios = append(ds.Trios, model.Trio{
			Child:  fields[0],
			Father: fields[1],
			Mother: fields[2],
		})
	}
	return scanner.Err()
}

// parsePhased reads the sectioned phased output file produced by output.WritePhased.
//
// Sections: #variants, #trios, #config, #phase, #blocks, #stats.
func parsePhased(path string) (*model.Dataset, *model.Result, config.Config, []string, error) {
	file, scanner, err := openLines(path)
	if err != nil {
		return nil, nil, config.Config{}, nil, err
	}
	defer file.Close()

	ds := &model.Dataset{GT: map[string][]model.Genotype{}}
	res := &model.Result{Samples: map[string]*model.SampleResult{}}
	cfg := config.Default()
	var notices []string
	section := ""
	headerSeen := map[string]bool{}

	flushSample := func(prev string, sites []model.SitePhase) {
		if prev == "" {
			return
		}
		sr := &model.SampleResult{Sample: prev, Sites: sites}
		res.Order = append(res.Order, prev)
		res.Samples[prev] = sr
	}

	var curSample string
	var curSites []model.SitePhase
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, "#") {
			flushSample(curSample, curSites)
			curSample, curSites = "", nil
			tokens := strings.Split(line, "\t")
			section = strings.TrimSpace(strings.TrimPrefix(tokens[0], "#"))
			continue
		}
		fields := strings.Split(line, "\t")
		if section != "config" && !headerSeen[section] {
			headerSeen[section] = true
			continue // column header row
		}
		switch section {
		case "variants":
			if len(fields) < 5 {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: malformed variants row", path, lineNo)
			}
			pos, convErr := strconv.Atoi(fields[2])
			if convErr != nil {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: invalid pos", path, lineNo)
			}
			ds.Variants = append(ds.Variants, model.Variant{
				ID: fields[0], Chrom: fields[1], Pos: pos, Ref: fields[3], Alt: fields[4],
			})
		case "trios":
			if len(fields) < 3 {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: malformed trios row", path, lineNo)
			}
			ds.Trios = append(ds.Trios, model.Trio{Child: fields[0], Father: fields[1], Mother: fields[2]})
		case "reads":
			if len(fields) < 4 {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: malformed reads row", path, lineNo)
			}
			if err := appendPhasedRead(ds, fields, path, lineNo); err != nil {
				return nil, nil, cfg, nil, err
			}
		case "config":
			if jsonErr := json.Unmarshal([]byte(line), &cfg); jsonErr != nil {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: malformed config: %w", path, lineNo, jsonErr)
			}
		case "phase":
			// sample variant_index variant_id status low_cov a0 a1 h1 h2 violation
			if len(fields) < 10 {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: malformed phase row", path, lineNo)
			}
			if _, ok := ds.GT[fields[0]]; !ok {
				ds.Samples = append(ds.Samples, fields[0])
				ds.GT[fields[0]] = make([]model.Genotype, len(ds.Variants))
			}
			if fields[0] != curSample {
				flushSample(curSample, curSites)
				curSample = fields[0]
				curSites = nil
			}
			sp, siteErr := parsePhaseRow(fields, path, lineNo)
			if siteErr != nil {
				return nil, nil, cfg, nil, siteErr
			}
			// Legacy 1.0 files recorded "." for both alleles of a violated
			// heterozygous site, which loses the genotype and cannot be
			// re-read correctly. Give a clear, locatable compatibility error
			// instead of an unlocatable verification failure.
			if sp.Violation && (sp.H1 == model.AlleleMissing || sp.H2 == model.AlleleMissing) {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: legacy phased file (format 1.0): violated site %q has missing allele calls; re-run phasing with the current version to produce a re-readable file", path, lineNo, sp.VariantID)
			}
			curSites = append(curSites, sp)
			ds.GT[fields[0]][sp.Index] = reconstructGenotype(sp)
		case "blocks":
			// sample block_id start end sites hets support conflict confidence error_rate status
			if len(fields) < 11 {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: malformed blocks row", path, lineNo)
			}
			b, blockErr := parseBlockRow(fields, path, lineNo)
			if blockErr != nil {
				return nil, nil, cfg, nil, blockErr
			}
			if sr, ok := res.Samples[fields[0]]; ok {
				sr.Blocks = append(sr.Blocks, b)
			}
		case "stats":
			if len(fields) < 2 {
				continue
			}
			if err := parseStatsRow(&res.Stats, fields[0], fields[1]); err != nil {
				return nil, nil, cfg, nil, fmt.Errorf("%s:%d: %w", path, lineNo, err)
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, nil, cfg, nil, err
	}
	flushSample(curSample, curSites)

	// Recompute per-sample counters from rows so reload yields an equal result.
	for _, sampleID := range res.Order {
		sr := res.Samples[sampleID]
		for _, sp := range sr.Sites {
			if sp.Hetero {
				sr.HetsTotal++
				if sp.Phased {
					sr.PhasedHets++
				} else {
					sr.UnphasedHets++
				}
				if sp.Violation {
					sr.ViolationHets++
				}
			}
		}
	}
	for i := range ds.Reads {
		sort.SliceStable(ds.Reads[i].Observations, func(a, b int) bool {
			return ds.Reads[i].Observations[a].VariantIdx < ds.Reads[i].Observations[b].VariantIdx
		})
	}
	return ds, res, cfg, notices, nil
}

// appendPhasedRead adds one #reads row, grouping rows into ds.Reads in
// deterministic (read_id, sample) first-appearance order.
func appendPhasedRead(ds *model.Dataset, f []string, path string, lineNo int) error {
	readID, sampleID, variantID, allele := f[0], f[1], f[2], f[3]
	idx := -1
	for i, variant := range ds.Variants {
		if variant.ID == variantID {
			idx = i
			break
		}
	}
	if idx < 0 {
		return fmt.Errorf("%s:%d: unknown variant id %q in reads", path, lineNo, variantID)
	}
	variant := ds.Variants[idx]
	if allele != variant.Ref && allele != variant.Alt {
		return fmt.Errorf("%s:%d: allele %q is not ref/alt for %s", path, lineNo, allele, variantID)
	}
	pos := -1
	for i := range ds.Reads {
		if ds.Reads[i].ID == readID && ds.Reads[i].Sample == sampleID {
			pos = i
			break
		}
	}
	if pos < 0 {
		ds.Reads = append(ds.Reads, model.Read{ID: readID, Sample: sampleID})
		pos = len(ds.Reads) - 1
	}
	ds.Reads[pos].Observations = append(ds.Reads[pos].Observations, model.ReadObservation{
		VariantIdx: idx, Allele: allele,
	})
	return nil
}

// parsePhaseRow parses one #phase row into a SitePhase.
func parsePhaseRow(f []string, path string, lineNo int) (model.SitePhase, error) {
	idx, err := strconv.Atoi(f[1])
	if err != nil {
		return model.SitePhase{}, fmt.Errorf("%s:%d: invalid variant index", path, lineNo)
	}
	sp := model.SitePhase{
		Index:     idx,
		VariantID: f[2],
		Phased:    f[3] == model.SiteStatusPhased,
		Hetero:    f[3] == model.SiteStatusPhased || f[3] == model.SiteStatusUnphased,
		LowCov:    f[4] == "1",
		A0:        f[5],
		A1:        f[6],
		Violation: f[9] == "1",
	}
	if f[7] == "A0" {
		sp.H1 = model.AlleleA0
	} else if f[7] == "A1" {
		sp.H1 = model.AlleleA1
	}
	if f[8] == "A0" {
		sp.H2 = model.AlleleA0
	} else if f[8] == "A1" {
		sp.H2 = model.AlleleA1
	}
	return sp, nil
}

// reconstructGenotype rebuilds an unphased Genotype from a phased site row.
func reconstructGenotype(sp model.SitePhase) model.Genotype {
	call := func(c model.AlleleCall) (string, bool) {
		switch c {
		case model.AlleleA0:
			return sp.A0, true
		case model.AlleleA1:
			return sp.A1, true
		default:
			return "", false
		}
	}
	a1, ok1 := call(sp.H1)
	a2, ok2 := call(sp.H2)
	switch {
	case !ok1 && !ok2:
		return model.Genotype{Present: true, Ploidy: 0}
	case !ok2:
		return model.Genotype{Present: true, Ploidy: 1, A1: a1}
	default:
		return model.Genotype{Present: true, Ploidy: 2, A1: a1, A2: a2}
	}
}

// parseBlockRow parses one #blocks row.
func parseBlockRow(f []string, path string, lineNo int) (model.PhaseBlock, error) {
	atoi := func(s string) (int, error) { return strconv.Atoi(s) }
	atof := func(s string) (float64, error) { return strconv.ParseFloat(s, 64) }
	start, err := atoi(f[2])
	if err != nil {
		return model.PhaseBlock{}, fmt.Errorf("%s:%d: invalid block start", path, lineNo)
	}
	end, err := atoi(f[3])
	if err != nil {
		return model.PhaseBlock{}, fmt.Errorf("%s:%d: invalid block end", path, lineNo)
	}
	sites, err := atoi(f[4])
	if err != nil {
		return model.PhaseBlock{}, err
	}
	hets, err := atoi(f[5])
	if err != nil {
		return model.PhaseBlock{}, err
	}
	support, err := atoi(f[6])
	if err != nil {
		return model.PhaseBlock{}, err
	}
	conflict, err := atoi(f[7])
	if err != nil {
		return model.PhaseBlock{}, err
	}
	confidence, err := atof(f[8])
	if err != nil {
		return model.PhaseBlock{}, err
	}
	errorRate, err := atof(f[9])
	if err != nil {
		return model.PhaseBlock{}, err
	}
	return model.PhaseBlock{
		ID: f[1], StartIdx: start, EndIdx: end, SiteCount: sites, HetsCount: hets,
		Support: support, Conflict: conflict, Confidence: confidence,
		ExpectedErrorRate: errorRate, Status: f[10],
	}, nil
}

// parseStatsRow maps a stats key to the RunStats field.
func parseStatsRow(st *model.RunStats, key, value string) error {
	atoi := func(s string) (int, error) { return strconv.Atoi(s) }
	var err error
	switch key {
	case "samples":
		st.Samples, err = atoi(value)
	case "variants":
		st.Variants, err = atoi(value)
	case "hets_total":
		st.HetsTotal, err = atoi(value)
	case "phased_hets":
		st.PhasedHets, err = atoi(value)
	case "unphased_hets":
		st.UnphasedHets, err = atoi(value)
	case "blocks_total":
		st.BlocksTotal, err = atoi(value)
	case "read_observations":
		st.ReadObservations, err = atoi(value)
	case "read_pairs":
		st.ReadPairs, err = atoi(value)
	case "read_edge_supports":
		st.ReadEdgeSupports, err = atoi(value)
	case "read_edge_conflicts":
		st.ReadEdgeConflicts, err = atoi(value)
	case "read_conflicts_adjudicated":
		st.ReadConflictsAdjudicated, err = atoi(value)
	case "trio_count":
		st.TrioCount, err = atoi(value)
	case "trio_constraint_sites":
		st.TrioConstraintSites, err = atoi(value)
	case "trio_violations":
		st.TrioViolations, err = atoi(value)
	case "breaks_by_threshold":
		st.BreaksByThreshold, err = atoi(value)
	case "low_coverage_sites":
		st.LowCoverageSites, err = atoi(value)
	case "missing_genotypes":
		st.MissingGenotypes, err = atoi(value)
	}
	return err
}

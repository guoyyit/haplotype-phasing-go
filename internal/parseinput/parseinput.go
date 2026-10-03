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
// The returned version is the format version recorded in the file header
// ("1.0" when absent); callers must re-emit the same version so that legacy
// result files verify byte-for-byte.
func LoadPhased(path string) (*model.Dataset, *model.Result, config.Config, []string, string, error) {
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
		gt, gtErr := parseGenotype(fields[2], ds.Variants[idx])
		if gtErr != nil {
			return fmt.Errorf("%s:%d: %w", path, lineNo, gtErr)
		}
		if !gt.Present {
			*notices = append(*notices, fmt.Sprintf("missing genotype: sample=%s variant=%s line=%d", sampleID, variantID, lineNo))
		}
		ds.GT[sampleID][idx] = gt
	}
	return scanner.Err()
}

// parseGenotype accepts ./. , ./A , A/. , A/A , A/B with either / or |
// separator. The missing token "." is symmetric on either allele position:
// one missing allele means a ploidy-1 partial genotype that keeps the single
// observed allele; two missing alleles mean a ploidy-0 full absence. An
// observed allele that is neither the site ref nor the site alt is rejected
// explicitly rather than silently carried into the output as a third allele.
func parseGenotype(token string, variant model.Variant) (model.Genotype, error) {
	token = strings.ReplaceAll(token, "|", "/")
	parts := strings.Split(token, "/")
	if len(parts) != 2 {
		return model.Genotype{}, fmt.Errorf("invalid genotype %q: expected allele/allele with / or | separator", token)
	}
	a1, a2 := strings.TrimSpace(parts[0]), strings.TrimSpace(parts[1])
	miss1, miss2 := a1 == "." || a1 == "", a2 == "." || a2 == ""
	if miss1 && miss2 {
		return model.Genotype{Present: true, Ploidy: 0}, nil
	}
	checkAllele := func(a string) error {
		if a != variant.Ref && a != variant.Alt {
			return fmt.Errorf("genotype allele %q at variant %s is neither ref %q nor alt %q",
				a, variant.ID, variant.Ref, variant.Alt)
		}
		return nil
	}
	switch {
	case miss1:
		if err := checkAllele(a2); err != nil {
			return model.Genotype{}, err
		}
		return model.Genotype{Present: true, Ploidy: 1, A1: a2}, nil
	case miss2:
		if err := checkAllele(a1); err != nil {
			return model.Genotype{}, err
		}
		return model.Genotype{Present: true, Ploidy: 1, A1: a1}, nil
	default:
		if err := checkAllele(a1); err != nil {
			return model.Genotype{}, err
		}
		if err := checkAllele(a2); err != nil {
			return model.Genotype{}, err
		}
		return model.Genotype{Present: true, Ploidy: 2, A1: a1, A2: a2}, nil
	}
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
func parsePhased(path string) (*model.Dataset, *model.Result, config.Config, []string, string, error) {
	file, scanner, err := openLines(path)
	if err != nil {
		return nil, nil, config.Config{}, nil, "", err
	}
	defer file.Close()

	ds := &model.Dataset{GT: map[string][]model.Genotype{}}
	res := &model.Result{Samples: map[string]*model.SampleResult{}}
	cfg := config.Default()
	var notices []string
	// readIndex maps (read_id, sample) to its position in ds.Reads so the
	// #reads section groups observations in O(1) instead of rescanning.
	readIndex := map[readKey]int{}
	variantIndex := map[string]int{}
	// Files without an explicit header are treated as the original 1.0
	// format; every released writer emits one.
	version := "1.0"
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
			if section == "haplophased" && len(tokens) >= 2 {
				version = strings.TrimSpace(tokens[1])
			}
			continue
		}
		fields := strings.Split(line, "\t")
		// The #stats section has no column-header row; every other
		// section skips its first non-comment row as its header.
		if section != "config" && section != "stats" && !headerSeen[section] {
			headerSeen[section] = true
			continue // column header row
		}
		switch section {
		case "variants":
			if len(fields) < 5 {
				return nil, nil, cfg, nil, version, fmt.Errorf("%s:%d: malformed variants row", path, lineNo)
			}
			pos, convErr := strconv.Atoi(fields[2])
			if convErr != nil {
				return nil, nil, cfg, nil, version, fmt.Errorf("%s:%d: invalid pos", path, lineNo)
			}
			ds.Variants = append(ds.Variants, model.Variant{
				ID: fields[0], Chrom: fields[1], Pos: pos, Ref: fields[3], Alt: fields[4],
			})
			variantIndex[fields[0]] = len(ds.Variants) - 1
		case "trios":
			if len(fields) < 3 {
				return nil, nil, cfg, nil, version, fmt.Errorf("%s:%d: malformed trios row", path, lineNo)
			}
			ds.Trios = append(ds.Trios, model.Trio{Child: fields[0], Father: fields[1], Mother: fields[2]})
		case "reads":
			if len(fields) < 4 {
				return nil, nil, cfg, nil, version, fmt.Errorf("%s:%d: malformed reads row", path, lineNo)
			}
			if err := appendPhasedRead(ds, readIndex, variantIndex, fields, path, lineNo); err != nil {
				return nil, nil, cfg, nil, version, err
			}
		case "config":
			if jsonErr := json.Unmarshal([]byte(line), &cfg); jsonErr != nil {
				return nil, nil, cfg, nil, version, fmt.Errorf("%s:%d: malformed config: %w", path, lineNo, jsonErr)
			}
		case "phase":
			// sample variant_index variant_id status low_cov a0 a1 h1 h2 violation
			if len(fields) < 10 {
				return nil, nil, cfg, nil, version, fmt.Errorf("%s:%d: malformed phase row", path, lineNo)
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
			sp, status, siteErr := parsePhaseRow(fields, path, lineNo)
			if siteErr != nil {
				return nil, nil, cfg, nil, version, siteErr
			}
			curSites = append(curSites, sp)
			ds.GT[fields[0]][sp.Index] = reconstructGenotype(sp, status, version == "1.0")
		case "blocks":
			// sample block_id start end sites hets support conflict confidence error_rate status
			if len(fields) < 11 {
				return nil, nil, cfg, nil, version, fmt.Errorf("%s:%d: malformed blocks row", path, lineNo)
			}
			b, blockErr := parseBlockRow(fields, path, lineNo)
			if blockErr != nil {
				return nil, nil, cfg, nil, version, blockErr
			}
			if sr, ok := res.Samples[fields[0]]; ok {
				sr.Blocks = append(sr.Blocks, b)
			}
		case "stats":
			if len(fields) < 2 {
				continue
			}
			if err := parseStatsRow(&res.Stats, fields[0], fields[1]); err != nil {
				return nil, nil, cfg, nil, version, fmt.Errorf("%s:%d: %w", path, lineNo, err)
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, nil, cfg, nil, version, err
	}
	flushSample(curSample, curSites)
	if version != "1.0" && version != "1.1" {
		return nil, nil, cfg, nil, version, fmt.Errorf("%s: unsupported phased format version %q", path, version)
	}

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
	return ds, res, cfg, notices, version, nil
}

// appendPhasedRead adds one #reads row, grouping rows into ds.Reads in
// deterministic (read_id, sample) first-appearance order.
func appendPhasedRead(ds *model.Dataset, readIndex map[readKey]int, variantIndex map[string]int, f []string, path string, lineNo int) error {
	readID, sampleID, variantID, allele := f[0], f[1], f[2], f[3]
	idx, known := variantIndex[variantID]
	if !known {
		return fmt.Errorf("%s:%d: unknown variant id %q in reads", path, lineNo, variantID)
	}
	variant := ds.Variants[idx]
	if allele != variant.Ref && allele != variant.Alt {
		return fmt.Errorf("%s:%d: allele %q is not ref/alt for %s", path, lineNo, allele, variantID)
	}
	key := readKey{sample: sampleID, readID: readID}
	pos, ok := readIndex[key]
	if !ok {
		ds.Reads = append(ds.Reads, model.Read{ID: readID, Sample: sampleID})
		pos = len(ds.Reads) - 1
		readIndex[key] = pos
	}
	ds.Reads[pos].Observations = append(ds.Reads[pos].Observations, model.ReadObservation{
		VariantIdx: idx, Allele: allele,
	})
	return nil
}

// parsePhaseRow parses one #phase row into a SitePhase.
func parsePhaseRow(f []string, path string, lineNo int) (model.SitePhase, string, error) {
	idx, err := strconv.Atoi(f[1])
	if err != nil {
		return model.SitePhase{}, "", fmt.Errorf("%s:%d: invalid variant index", path, lineNo)
	}
	status := f[3]
	sp := model.SitePhase{
		Index:     idx,
		VariantID: f[2],
		Phased:    status == model.SiteStatusPhased,
		// UNPHASED covers both violating heterozygotes and violating
		// homozygotes; the latter are distinguished below by their two
		// identical haplotype calls (or, in 1.0 files, by reconstruction).
		Hetero: status == model.SiteStatusPhased ||
			(status == model.SiteStatusUnphased && f[7] != f[8]),
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
	return sp, status, nil
}

// reconstructGenotype rebuilds an unphased Genotype from a phased site row.
// legacyV1Dot0 records whether the row came from a format 1.0 file, where
// violating heterozygous sites were written with two "." calls. In that
// case a dots-only UNPHASED row is reconstructed as a heterozygous genotype,
// which is what re-deriving the violation and re-emitting the row verbatim
// requires; format 1.1 distinguishes hetero (A0/A1) from homozygous
// violations explicitly.
func reconstructGenotype(sp model.SitePhase, status string, legacyV1Dot0 bool) model.Genotype {
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
		if legacyV1Dot0 && status == model.SiteStatusUnphased && sp.Violation {
			// Format 1.0 encoded violating heterozygous sites as two
			// missing calls; reconstruct the hetero genotype so the
			// violation re-derives identically and the 1.0 row re-emits.
			return model.Genotype{Present: true, Ploidy: 2, A1: sp.A0, A2: sp.A1}
		}
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

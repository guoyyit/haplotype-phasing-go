package output

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"strconv"

	"github.com/guoyyit/haplotype-phasing-go/internal/config"
	"github.com/guoyyit/haplotype-phasing-go/internal/model"
)

// FileVersion is the stable phased-output format version.
const FileVersion = "1.0"

// WritePhased writes the stable, re-readable phased result format.
func WritePhased(w io.Writer, ds *model.Dataset, res *model.Result, cfg config.Config) error {
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	fmt.Fprintf(bw, "#haplophased\t%s\n", FileVersion)

	// #variants
	fmt.Fprintln(bw, "#variants")
	fmt.Fprintln(bw, "variant_id\tchrom\tpos\tref\talt")
	for _, v := range ds.Variants {
		fmt.Fprintf(bw, "%s\t%s\t%d\t%s\t%s\n", v.ID, v.Chrom, v.Pos, v.Ref, v.Alt)
	}

	// #trios
	fmt.Fprintln(bw, "#trios")
	fmt.Fprintln(bw, "child\tfather\tmother")
	for _, t := range ds.Trios {
		fmt.Fprintf(bw, "%s\t%s\t%s\n", t.Child, t.Father, t.Mother)
	}

	// #reads (required to reproduce read-pair evidence on re-read)
	fmt.Fprintln(bw, "#reads")
	fmt.Fprintln(bw, "read_id\tsample\tvariant_id\tallele")
	for _, read := range ds.Reads {
		for _, obs := range read.Observations {
			fmt.Fprintf(bw, "%s\t%s\t%s\t%s\n", read.ID, read.Sample, ds.Variants[obs.VariantIdx].ID, obs.Allele)
		}
	}

	// #config (single compact JSON line)
	fmt.Fprintln(bw, "#config")
	raw, err := json.Marshal(cfg)
	if err != nil {
		return err
	}
	fmt.Fprintln(bw, string(raw))

	// #phase
	fmt.Fprintln(bw, "#phase")
	fmt.Fprintln(bw, "sample\tvariant_index\tvariant_id\tstatus\tlow_cov\ta0\ta1\th1\th2\tviolation")
	for _, sample := range res.Order {
		sr := res.Samples[sample]
		for _, sp := range sr.Sites {
			fmt.Fprintf(bw, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				sample, sp.Index, sp.VariantID, siteStatus(sp),
				bit(sp.LowCov), sp.A0, sp.A1,
				alleleName(sp.H1, sp.A0, sp.A1),
				alleleName(sp.H2, sp.A0, sp.A1),
				bit(sp.Violation))
		}
	}

	// #blocks
	fmt.Fprintln(bw, "#blocks")
	fmt.Fprintln(bw, "sample\tblock_id\tstart_idx\tend_idx\tsite_count\thets_count\tsupport\tconflict\tconfidence\terror_rate\tstatus")
	for _, sample := range res.Order {
		sr := res.Samples[sample]
		for _, b := range sr.Blocks {
			fmt.Fprintf(bw, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%s\n",
				sample, b.ID, b.StartIdx, b.EndIdx, b.SiteCount, b.HetsCount,
				b.Support, b.Conflict, b.Confidence, b.ExpectedErrorRate, b.Status)
		}
	}

	// #stats
	fmt.Fprintln(bw, "#stats")
	writeStats(bw, res.Stats)
	return bw.Flush()
}

// siteStatus maps a SitePhase to the documented status token.
func siteStatus(sp model.SitePhase) string {
	switch {
	case sp.Violation:
		return model.SiteStatusUnphased
	case sp.Phased:
		return model.SiteStatusPhased
	case sp.Hetero:
		return model.SiteStatusUnphased
	default:
		if sp.H1 == model.AlleleMissing {
			return model.SiteStatusMissing
		}
		if sp.H2 == model.AlleleMissing {
			return model.SiteStatusPartial
		}
		return model.SiteStatusHomozygous
	}
}

func alleleName(call model.AlleleCall, a0, a1 string) string {
	switch call {
	case model.AlleleA0:
		return "A0"
	case model.AlleleA1:
		return "A1"
	default:
		return "."
	}
}

func bit(v bool) string {
	if v {
		return "1"
	}
	return "0"
}

func writeStats(bw *bufio.Writer, st model.RunStats) {
	rows := [][2]string{
		{"samples", strconv.Itoa(st.Samples)},
		{"variants", strconv.Itoa(st.Variants)},
		{"hets_total", strconv.Itoa(st.HetsTotal)},
		{"phased_hets", strconv.Itoa(st.PhasedHets)},
		{"unphased_hets", strconv.Itoa(st.UnphasedHets)},
		{"blocks_total", strconv.Itoa(st.BlocksTotal)},
		{"read_observations", strconv.Itoa(st.ReadObservations)},
		{"read_pairs", strconv.Itoa(st.ReadPairs)},
		{"read_edge_supports", strconv.Itoa(st.ReadEdgeSupports)},
		{"read_edge_conflicts", strconv.Itoa(st.ReadEdgeConflicts)},
		{"read_conflicts_adjudicated", strconv.Itoa(st.ReadConflictsAdjudicated)},
		{"trio_count", strconv.Itoa(st.TrioCount)},
		{"trio_constraint_sites", strconv.Itoa(st.TrioConstraintSites)},
		{"trio_violations", strconv.Itoa(st.TrioViolations)},
		{"breaks_by_threshold", strconv.Itoa(st.BreaksByThreshold)},
		{"low_coverage_sites", strconv.Itoa(st.LowCoverageSites)},
		{"missing_genotypes", strconv.Itoa(st.MissingGenotypes)},
	}
	for _, row := range rows {
		fmt.Fprintf(bw, "%s\t%s\n", row[0], row[1])
	}
}

// WriteReport writes a human-readable summary report.
func WriteReport(w io.Writer, ds *model.Dataset, res *model.Result, cfg config.Config) error {
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	fmt.Fprintln(bw, "Haplotype phasing report")
	fmt.Fprintf(bw, "samples=%d variants=%d\n", res.Stats.Samples, res.Stats.Variants)
	fmt.Fprintf(bw, "hets_total=%d phased_hets=%d unphased_hets=%d blocks=%d\n",
		res.Stats.HetsTotal, res.Stats.PhasedHets, res.Stats.UnphasedHets, res.Stats.BlocksTotal)
	fmt.Fprintf(bw, "read_observations=%d read_pairs=%d read_edge_conflicts=%d adjudicated=%d\n",
		res.Stats.ReadObservations, res.Stats.ReadPairs, res.Stats.ReadEdgeConflicts, res.Stats.ReadConflictsAdjudicated)
	fmt.Fprintf(bw, "trios=%d trio_constraint_sites=%d trio_violations=%d\n",
		res.Stats.TrioCount, res.Stats.TrioConstraintSites, res.Stats.TrioViolations)
	fmt.Fprintf(bw, "breaks_by_threshold=%d low_coverage_sites=%d missing_genotypes=%d\n",
		res.Stats.BreaksByThreshold, res.Stats.LowCoverageSites, res.Stats.MissingGenotypes)
	fmt.Fprintln(bw)
	for _, sample := range res.Order {
		sr := res.Samples[sample]
		fmt.Fprintf(bw, "sample %s: hets=%d phased=%d unphased=%d violations=%d blocks=%d\n",
			sample, sr.HetsTotal, sr.PhasedHets, sr.UnphasedHets, sr.ViolationHets, len(sr.Blocks))
		for _, b := range sr.Blocks {
			fmt.Fprintf(bw, "  %s sites=%d-%d hets=%d support=%d conflict=%d confidence=%.4f error=%.4f %s\n",
				b.ID, b.StartIdx, b.EndIdx, b.HetsCount, b.Support, b.Conflict,
				b.Confidence, b.ExpectedErrorRate, b.Status)
		}
	}
	return bw.Flush()
}

package phasing

// Variant describes one biallelic diploid autosomal variant site.
type Variant struct {
	Index int    `json:"index"`
	ID    string `json:"id"`
	Chrom string `json:"chrom"`
	Pos   int64  `json:"pos"`
	Ref   string `json:"ref"`
	Alt   string `json:"alt"`
}

// GenotypeCall is an unordered diploid genotype. Allele values are 0 (ref),
// 1 (alt) or -1 when the allele is missing (partial deletion, e.g. "./1").
type GenotypeCall struct {
	Present bool `json:"-"`
	A1      int8 `json:"-"`
	A2      int8 `json:"-"`
}

// Sample holds the per-variant genotype calls of one sample in variant order.
type Sample struct {
	ID string         `json:"id"`
	GT []GenotypeCall `json:"-"`
}

// SiteObs is one allele observation of a sequencing read at one variant site.
type SiteObs struct {
	Variant int     `json:"variant_idx"`
	Allele  int8    `json:"allele"`
	Weight  float64 `json:"weight"`
}

// Read groups all observations belonging to the same physical sequencing read.
type Read struct {
	ID     string    `json:"id"`
	Sample string    `json:"sample_id"`
	Obs    []SiteObs `json:"obs"`
}

// Trio is an optional child/father/mother pedigree annotation.
type Trio struct {
	Child  string `json:"child"`
	Father string `json:"father"`
	Mother string `json:"mother"`
}

// Input bundles everything the phasing engine consumes.
type Input struct {
	Variants []Variant
	Samples  []Sample
	Reads    []Read
	Trio     *Trio
	Config   Config
}

// SiteStatus values mirror the documented output contract.
const (
	SitePhased   = "phased"
	SiteUnphased = "unphased"
	SiteMissing  = "missing"
)

// SitePhase is the per-sample, per-site result entry.
type SitePhase struct {
	VariantIdx      int    `json:"variant_idx"`
	VariantID       string `json:"variant_id"`
	Status          string `json:"status"`
	H1              int8   `json:"h1"`
	H2              int8   `json:"h2"`
	BlockID         int    `json:"block_id"`
	PaternalH1      bool   `json:"paternal_h1"`
	OriginResolved  bool   `json:"origin_resolved"`
	MendelViolation bool   `json:"mendel_violation"`
	LowCoverage     bool   `json:"low_coverage"`
	ReadSupport     int    `json:"read_support"`
}

// BlockInfo describes one phased block of one sample.
type BlockInfo struct {
	SampleID          string  `json:"sample_id"`
	BlockID           int     `json:"block_id"`
	StartVariantIdx   int     `json:"start_variant_idx"`
	EndVariantIdx     int     `json:"end_variant_idx"`
	Hets              int     `json:"heterozygous_sites"`
	SupportingReads   int     `json:"supporting_reads"`
	Confidence        float64 `json:"confidence"`
	ExpectedErrorRate float64 `json:"expected_error_rate"`
}

// ConflictEntry records one adjudicated conflicting observation.
type ConflictEntry struct {
	SampleID    string `json:"sample_id"`
	ReadID      string `json:"read_id,omitempty"`
	VariantI    int    `json:"variant_i"`
	VariantJ    int    `json:"variant_j"`
	VotedParity int    `json:"voted_parity"`
	KeptParity  int    `json:"kept_parity"`
	Kind        string `json:"kind"`
}

// ConflictReport aggregates all evidence conflicts.
type ConflictReport struct {
	ReadReadConflicts     int             `json:"read_read_conflicts"`
	ReadTrioConflicts     int             `json:"read_trio_conflicts"`
	ReadGenotypeConflicts int             `json:"read_genotype_conflicts"`
	TrioMendelViolations  int             `json:"trio_mendel_violations"`
	Entries               []ConflictEntry `json:"entries"`
}

// SampleResult is the phasing output for one sample.
type SampleResult struct {
	ID     string      `json:"sample_id"`
	Blocks []BlockInfo `json:"blocks"`
	Sites  []SitePhase `json:"sites"`
}

// Result is the complete, serializable phasing result.
type Result struct {
	SchemaVersion string         `json:"schema_version"`
	Config        ResolvedConfig `json:"config"`
	Variants      []Variant      `json:"variants"`
	Samples       []SampleResult `json:"samples"`
	Conflicts     ConflictReport `json:"conflicts"`
	Stats         RunStats       `json:"stats"`
}

// RunStats summarizes the run.
type RunStats struct {
	Samples           int `json:"samples"`
	VariantSites      int `json:"variant_sites"`
	HeterozygousSites int `json:"heterozygous_sites"`
	PhasedHets        int `json:"phased_hets"`
	UnphasedHets      int `json:"unphased_hets"`
	PartialSites      int `json:"partial_sites"`
	MissingSites      int `json:"missing_sites"`
	LowCoverageSites  int `json:"low_coverage_sites"`
	CorrectlyPhased   int `json:"correctly_phased"`
}

package model

// AlleleCall identifies an allele relative to a site's canonical allele order.
type AlleleCall int

const (
	AlleleMissing AlleleCall = 0
	AlleleA0      AlleleCall = 1
	AlleleA1      AlleleCall = 2
)

// Stable status tokens used in the phased output file.
const (
	SiteStatusPhased     = "PHASED"
	SiteStatusUnphased   = "UNPHASED"
	SiteStatusHomozygous = "HOMOZYGOUS"
	SiteStatusPartial    = "PARTIAL"
	SiteStatusMissing    = "MISSING"
	BlockStatusOK        = "OK"
	BlockStatusLowConf   = "LOW_CONF"
	BlockStatusSingleton = "SINGLETON"
)

// Variant describes one diploid autosomal biallelic variant site.
type Variant struct {
	ID    string
	Chrom string
	Pos   int
	Ref   string
	Alt   string
}

// Genotype is an unphased, possibly partially missing diploid genotype.
// Ploidy is the number of observed alleles (0, 1 or 2).
type Genotype struct {
	Present bool
	Ploidy  int
	A1      string
	A2      string
}

// ReadObservation is one allele observation of a read at one site.
type ReadObservation struct {
	VariantIdx int
	Allele     string
}

// Read is a sequencing read spanning one or more sites of one sample.
type Read struct {
	ID           string
	Sample       string
	Observations []ReadObservation
}

// Trio is a father-mother-child family.
type Trio struct {
	Child  string
	Father string
	Mother string
}

// Dataset is the complete parsed input.
type Dataset struct {
	Variants []Variant
	Samples  []string
	GT       map[string][]Genotype
	Reads    []Read
	Trios    []Trio
}

// SitePhase is the phasing call for one sample at one site.
type SitePhase struct {
	Index     int
	VariantID string
	Phased    bool
	Hetero    bool
	LowCov    bool
	A0        string
	A1        string
	H1        AlleleCall
	H2        AlleleCall
	Violation bool
}

// PhaseBlock is a maximal set of adjacent sites connected by sufficient evidence.
type PhaseBlock struct {
	ID                string
	StartIdx          int
	EndIdx            int
	SiteCount         int
	HetsCount         int
	Support           int
	Conflict          int
	Confidence        float64
	ExpectedErrorRate float64
	Status            string
}

// SampleResult is the full per-sample phasing output.
type SampleResult struct {
	Sample        string
	Sites         []SitePhase
	Blocks        []PhaseBlock
	HetsTotal     int
	PhasedHets    int
	UnphasedHets  int
	ViolationHets int
}

// RunStats summarizes a phasing run.
type RunStats struct {
	Samples                  int
	Variants                 int
	HetsTotal                int
	PhasedHets               int
	UnphasedHets             int
	BlocksTotal              int
	ReadObservations         int
	ReadPairs                int
	ReadEdgeSupports         int
	ReadEdgeConflicts        int
	ReadConflictsAdjudicated int
	TrioCount                int
	TrioConstraintSites      int
	TrioViolations           int
	BreaksByThreshold        int
	LowCoverageSites         int
	MissingGenotypes         int
}

// Result is a complete phasing result set.
type Result struct {
	Samples map[string]*SampleResult
	Order   []string
	Stats   RunStats
}

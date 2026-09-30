package phasing

import (
	"errors"
	"fmt"
)

// Config declares every user-tunable option (see docs for semantics).
//
// Block breaking rule: a read-derived adjacency between two consecutive
// heterozygous sites is accepted only when the number of adjudicated
// supporting votes is >= MinEdgeSupport AND the posterior edge confidence is
// >= MinBlockConfidence. Every rejected adjacency deterministically breaks the
// phasing block.
type Config struct {
	// MinEdgeSupport is the minimum number of independent read votes (after
	// conflict adjudication) required to join two sites into one block.
	MinEdgeSupport int `json:"min_edge_support"`

	// MinReadsPerSite is the coverage floor: sites covered by fewer reads are
	// still phased when valid evidence exists but are flagged low_coverage.
	MinReadsPerSite int `json:"min_reads_per_site"`

	// MinBlockConfidence is the minimum posterior confidence (in [0.5, 1])
	// required to accept a read-derived adjacency and reported per block.
	MinBlockConfidence float64 `json:"min_block_confidence"`

	// ReadErrorRate is the per-observation sequencing/allele error probability
	// used for posterior confidence calculations and conflict adjudication.
	ReadErrorRate float64 `json:"read_error_rate"`

	// SeedEnabled controls tie resolution. When true, parity votes that tie
	// are resolved with the deterministic PRNG seeded by RandomSeed. When
	// false, ties are resolved by the lexicographically smallest read ID.
	SeedEnabled bool  `json:"seed_enabled"`
	RandomSeed  int64 `json:"random_seed"`
}

// Default configuration values, documented and stable.
const (
	defaultMinEdgeSupport     = 1
	defaultMinReadsPerSite    = 1
	defaultMinBlockConfidence = 0.90
	defaultReadErrorRate      = 0.01
)

// DefaultConfig returns the documented default configuration.
func DefaultConfig() Config {
	return Config{
		MinEdgeSupport:     defaultMinEdgeSupport,
		MinReadsPerSite:    defaultMinReadsPerSite,
		MinBlockConfidence: defaultMinBlockConfidence,
		ReadErrorRate:      defaultReadErrorRate,
		SeedEnabled:        false,
		RandomSeed:         0,
	}
}

// ResolvedConfig is the fully defaulted, validated configuration stored in output.
type ResolvedConfig struct {
	Config
}

// Resolve fills defaults for unset fields and validates the configuration.
// A zero Config resolves to DefaultConfig.
func (c Config) Resolve() (ResolvedConfig, error) {
	d := DefaultConfig()
	if c.MinEdgeSupport == 0 {
		c.MinEdgeSupport = d.MinEdgeSupport
	}
	if c.MinReadsPerSite == 0 {
		c.MinReadsPerSite = d.MinReadsPerSite
	}
	if c.MinBlockConfidence == 0 {
		c.MinBlockConfidence = d.MinBlockConfidence
	}
	if c.ReadErrorRate == 0 {
		c.ReadErrorRate = d.ReadErrorRate
	}

	var errs []error
	if c.MinEdgeSupport < 1 {
		errs = append(errs, errors.New("min_edge_support must be >= 1"))
	}
	if c.MinReadsPerSite < 1 {
		errs = append(errs, errors.New("min_reads_per_site must be >= 1"))
	}
	if c.ReadErrorRate <= 0 || c.ReadErrorRate >= 0.5 {
		errs = append(errs, errors.New("read_error_rate must be in the open interval (0, 0.5)"))
	}
	if c.MinBlockConfidence < 0.5 || c.MinBlockConfidence > 1.0 {
		errs = append(errs, errors.New("min_block_confidence must be in [0.5, 1.0]"))
	}
	if len(errs) > 0 {
		return ResolvedConfig{}, fmt.Errorf("invalid configuration: %w", errs[0])
	}
	return ResolvedConfig{Config: c}, nil
}

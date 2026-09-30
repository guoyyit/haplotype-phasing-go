package config

import (
	"encoding/json"
	"fmt"
	"os"
)

// Config holds all user-tunable, documented phasing parameters.
type Config struct {
	// MinEdgeSupport is the minimum number of read-pair votes on an edge
	// between two adjacent heterozygous sites for the edge to connect them
	// into one phase block. 0 disables the depth requirement.
	MinEdgeSupport int `json:"min_edge_support"`
	// MinEdgeStrength is the minimum fraction of votes agreeing with the
	// majority orientation. An edge whose strength is below this value
	// breaks the phase block.
	MinEdgeStrength float64 `json:"min_edge_strength"`
	// MinBlockConfidence flags a block as LOW_CONF when its confidence is
	// below this value; the block is still reported.
	MinBlockConfidence float64 `json:"min_block_confidence"`
	// BlockMaxErrorRate is the maximum allowed per-edge expected error
	// rate; edges above it break the phase block.
	BlockMaxErrorRate float64 `json:"block_max_error_rate"`
	// MinSiteDepth marks a site LOW_COV when fewer reads cover it. It is
	// counted and reported but never aborts processing.
	MinSiteDepth int `json:"min_site_depth"`
	// Seed drives the deterministic tie-break of equally split evidence.
	Seed int64 `json:"seed"`
	// EnableTrio turns trio-inheritance evidence on/off.
	EnableTrio bool `json:"enable_trio"`
	// EnableReads turns read-pair evidence on/off.
	EnableReads bool `json:"enable_reads"`
}

// Default returns the default configuration.
func Default() Config {
	return Config{
		MinEdgeSupport:     2,
		MinEdgeStrength:    0.8,
		MinBlockConfidence: 0.95,
		BlockMaxErrorRate:  0.25,
		MinSiteDepth:       3,
		Seed:               1,
		EnableTrio:         true,
		EnableReads:        true,
	}
}

// Validate checks all configuration ranges.
func (config Config) Validate() error {
	if config.MinEdgeSupport < 0 {
		return fmt.Errorf("min_edge_support must be >= 0, got %d", config.MinEdgeSupport)
	}
	if config.MinEdgeStrength < 0.5 || config.MinEdgeStrength > 1 {
		return fmt.Errorf("min_edge_strength must be in [0.5,1], got %v", config.MinEdgeStrength)
	}
	if config.MinBlockConfidence <= 0 || config.MinBlockConfidence > 1 {
		return fmt.Errorf("min_block_confidence must be in (0,1], got %v", config.MinBlockConfidence)
	}
	if config.BlockMaxErrorRate < 0 || config.BlockMaxErrorRate > 1 {
		return fmt.Errorf("block_max_error_rate must be in [0,1], got %v", config.BlockMaxErrorRate)
	}
	if config.MinSiteDepth < 0 {
		return fmt.Errorf("min_site_depth must be >= 0, got %d", config.MinSiteDepth)
	}
	return nil
}

// Load reads a JSON configuration file; missing path or file means defaults.
func Load(path string) (config Config, err error) {
	config = Default()
	if path == "" {
		return config, nil
	}
	data, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return config, nil
		}
		return config, err
	}
	if err := json.Unmarshal(data, &config); err != nil {
		return config, fmt.Errorf("parsing config %s: %w", path, err)
	}
	if err := config.Validate(); err != nil {
		return config, err
	}
	return config, nil
}

package config

import (
	"os"
	"path/filepath"
	"testing"
)

func TestDefaultsValid(t *testing.T) {
	if err := Default().Validate(); err != nil {
		t.Fatalf("default config must be valid: %v", err)
	}
}

func TestLoadMissingFileUsesDefaults(t *testing.T) {
	cfg, err := Load(filepath.Join(t.TempDir(), "absent.json"))
	if err != nil {
		t.Fatal(err)
	}
	if cfg != Default() {
		t.Fatal("absent config file must fall back to defaults")
	}
}

func TestLoadOverrideAndReject(t *testing.T) {
	dir := t.TempDir()
	good := filepath.Join(dir, "good.json")
	if err := os.WriteFile(good, []byte(`{"min_edge_support": 4, "seed": 7}`), 0644); err != nil {
		t.Fatal(err)
	}
	cfg, err := Load(good)
	if err != nil {
		t.Fatal(err)
	}
	if cfg.MinEdgeSupport != 4 || cfg.Seed != 7 {
		t.Fatalf("override not applied: %+v", cfg)
	}

	bad := filepath.Join(dir, "bad.json")
	if err := os.WriteFile(bad, []byte(`{"min_edge_strength": 0.1}`), 0644); err != nil {
		t.Fatal(err)
	}
	if _, err := Load(bad); err == nil {
		t.Fatal("strength below 0.5 must be rejected")
	}
}

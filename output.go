package phasing

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"os"
)

// SchemaVersion identifies the stable output format.
const SchemaVersion = "haplotyper-result/v1"

// WriteResult writes a result to a local file in the stable output format.
func WriteResult(path string, r *Result) error {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetIndent("", "  ")
	enc.SetEscapeHTML(false)
	if err := enc.Encode(r); err != nil {
		return fmt.Errorf("encode result: %w", err)
	}
	// Trailing newline for friendly local-file usage; bytes are otherwise
	// identical run to run because every field is populated in fixed order.
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		return fmt.Errorf("write result: %w", err)
	}
	return nil
}

// ReadResult reads a result previously produced by WriteResult.
func ReadResult(path string) (*Result, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open result: %w", err)
	}
	defer f.Close()
	return DecodeResult(f)
}

// ReadResultFile returns the raw bytes of a written result file.
func ReadResultFile(path string) ([]byte, error) {
	return os.ReadFile(path)
}

// DecodeResult parses a result from a reader in the stable output format.
func DecodeResult(r io.Reader) (*Result, error) {
	dec := json.NewDecoder(r)
	var res Result
	if err := dec.Decode(&res); err != nil {
		return nil, fmt.Errorf("decode result: %w", err)
	}
	if res.SchemaVersion != SchemaVersion {
		return nil, fmt.Errorf("unsupported schema_version %q, want %q", res.SchemaVersion, SchemaVersion)
	}
	return &res, nil
}

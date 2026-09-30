package phasing

import (
	"math"
	"testing"
)

func TestEdgeConfidence(t *testing.T) {
	e := 0.01
	// Unanimous single read: posterior of winning parity is 1-e.
	if got := edgeConfidence(1, 0, e); math.Abs(got-0.99) > 1e-9 {
		t.Fatalf("unanimous confidence = %v, want 0.99", got)
	}
	// Balanced evidence: exactly 0.5 confidence for either parity.
	if got := edgeConfidence(2, 2, e); math.Abs(got-0.5) > 1e-12 {
		t.Fatalf("tie confidence = %v, want 0.5", got)
	}
	// Strong disagreement resolves decisively.
	if got := edgeConfidence(5, 0, e); !(got > 0.99999) {
		t.Fatalf("5 unanimous votes confidence = %v, want ~1", got)
	}
	// No votes carries no confidence.
	if got := edgeConfidence(0, 0, e); got != 0 {
		t.Fatalf("empty edge confidence = %v, want 0", got)
	}
}

func TestRNGDeterministic(t *testing.T) {
	a := newRNG(12345)
	b := newRNG(12345)
	for i := 0; i < 1000; i++ {
		if a.next64() != b.next64() {
			t.Fatal("PCG streams diverged for the same seed")
		}
	}
	// intn stays in range.
	r := newRNG(7)
	for i := 0; i < 1000; i++ {
		if n := r.intn(3); n < 0 || n >= 3 {
			t.Fatalf("intn out of range: %d", n)
		}
	}
}

func TestConfigRejectsInvalid(t *testing.T) {
	bad := []Config{
		{MinEdgeSupport: -1},
		{ReadErrorRate: 0.6},
		{MinBlockConfidence: 0.4},
		{MinReadsPerSite: -1},
	}
	for i, c := range bad {
		if _, err := c.Resolve(); err == nil {
			t.Fatalf("bad config #%d accepted: %+v", i, c)
		}
	}
	rc, err := DefaultConfig().Resolve()
	if err != nil {
		t.Fatalf("default config invalid: %v", err)
	}
	if rc.MinBlockConfidence != 0.90 || rc.ReadErrorRate != 0.01 {
		t.Fatalf("unexpected defaults: %+v", rc)
	}
}

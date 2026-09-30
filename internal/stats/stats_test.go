package stats

import (
	"math"
	"testing"
)

func TestEdgeErrorRate(t *testing.T) {
	if got := EdgeErrorRate(0, 0); math.Abs(got-0.5) > 1e-9 {
		t.Fatalf("no evidence should be 0.5, got %v", got)
	}
	// Two unanimous votes: Jeffreys smoothed (0.5)/(3)=0.1667.
	if got := EdgeErrorRate(2, 0); math.Abs(got-1.0/6.0) > 1e-9 {
		t.Fatalf("expected 1/6, got %v", got)
	}
	if EdgeErrorRate(3, 1) >= EdgeErrorRate(1, 3) {
		t.Fatal("higher conflict fraction must raise error rate")
	}
}

func TestBlockErrorRate(t *testing.T) {
	if got := BlockErrorRate(nil); got != 0 {
		t.Fatalf("empty block error must be 0, got %v", got)
	}
	e := 0.1
	want := 1 - (1-e)*(1-e)
	if got := BlockErrorRate([]float64{e, e}); math.Abs(got-want) > 1e-12 {
		t.Fatalf("expected %v, got %v", want, got)
	}
	if got := BlockConfidence([]float64{0, 0}); got != 1 {
		t.Fatalf("perfect edges must give confidence 1, got %v", got)
	}
}

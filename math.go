package phasing

import "math"

// edgeConfidence returns the posterior probability that the adjacency parity
// is "flip" given nFlip flip votes and nSame same-orientation votes, under a
// fixed per-observation error rate e.
//
// A vote is correct with probability 1-e and flipped with probability e, with
// a uniform prior over the two parities, so:
//
//	P(flip|votes) proportional to (1-e)^nFlip e^nSame
//	P(same|votes) proportional to (1-e)^nSame e^Flip
//
// Computed in log space for numerical stability; the returned confidence is
// always the probability of the winning parity (>= 0.5).
func edgeConfidence(nFlip, nSame int, e float64) float64 {
	if nFlip == 0 && nSame == 0 {
		return 0
	}
	logFlip := float64(nFlip)*math.Log1p(-e) + float64(nSame)*math.Log(e)
	logSame := float64(nSame)*math.Log1p(-e) + float64(nFlip)*math.Log(e)
	// Normalize with the max log term to avoid underflow.
	m := math.Max(logFlip, logSame)
	pFlip := math.Exp(logFlip - m)
	pSame := math.Exp(logSame - m)
	p := pFlip / (pFlip + pSame)
	if p < 0.5 {
		p = 1 - p
	}
	return p
}

// rng64 is a deterministic permuted congruential generator (PCG-XSH-RR 64/32
// state, emitting 64 bits from two 32-bit draws). It exists so that seed-based
// tie resolution is fully reproducible and independent of the host's rand
// version.
type rng64 struct {
	state uint64
	inc   uint64
}

func newRNG(seed int64) *rng64 {
	r := &rng64{state: 0, inc: 0x853c49e6748fea9b}
	r.next32()
	r.state += uint64(seed)
	r.next32()
	return r
}

const pcgMultiplier = 6364136223846793005

func (r *rng64) next32() uint32 {
	old := r.state
	r.state = old*pcgMultiplier + r.inc
	xorshifted := uint32(((old >> 18) ^ old) >> 27)
	rot := uint32(old >> 59)
	return (xorshifted >> rot) | (xorshifted << ((-rot) & 31))
}

func (r *rng64) next64() uint64 {
	return uint64(r.next32())<<32 | uint64(r.next32())
}

// intn returns a deterministic value in [0, n).
func (r *rng64) intn(n int) int {
	if n <= 0 {
		return 0
	}
	return int(r.next64() % uint64(n))
}

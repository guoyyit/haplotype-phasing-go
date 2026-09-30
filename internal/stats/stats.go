package stats

// priorAlpha is the symmetric Beta(alpha,alpha) prior pseudo-count used
// for the Laplace/Jeffreys-smoothed edge error estimate.
const priorAlpha = 0.5

// EdgeErrorRate computes the per-edge expected switching-error probability.
//
// Semantics: each read-pair vote is a Bernoulli trial of the majority
// orientation. With a symmetric Beta(0.5,0.5) prior the posterior mean of
// the probability of a disagreeing (switching) observation is
//
//	e = (conflict + 0.5) / (support + conflict + 1)
//
// With zero votes the edge is maximally uncertain (e = 0.5).
func EdgeErrorRate(support, conflict int) float64 {
	if support < 0 {
		support = 0
	}
	if conflict < 0 {
		conflict = 0
	}
	return (float64(conflict) + priorAlpha) / (float64(support+conflict) + 2*priorAlpha)
}

// BlockErrorRate aggregates edge error probabilities over a block.
//
// Edges are treated as independent orientation decisions, so the block's
// expected switching-event rate is 1 - product(1 - e_i): the probability
// that at least one internal edge is flipped. An empty list yields 0.
func BlockErrorRate(edgeRates []float64) float64 {
	product := 1.0
	for _, rate := range edgeRates {
		if rate < 0 {
			rate = 0
		}
		if rate > 1 {
			rate = 1
		}
		product *= 1 - rate
	}
	rate := 1 - product
	if rate < 0 {
		return 0
	}
	if rate > 1 {
		return 1
	}
	return rate
}

// BlockConfidence returns 1 - BlockErrorRate.
func BlockConfidence(edgeRates []float64) float64 {
	return 1 - BlockErrorRate(edgeRates)
}

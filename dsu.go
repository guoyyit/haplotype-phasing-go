package phasing

// DSU is a weighted (parity) disjoint-set union over phase variables.
// Each variable x carries a bit phase[x]. A recorded constraint has the form
// phase[x] XOR phase[y] = p (p=0: same orientation, p=1: flipped orientation).
// parity[x] stores phase[x] XOR phase[parent(x)].
type DSU struct {
	parent []int
	parity []int8
	rank   []int
}

// NewDSU creates an isolated parity-DSU over n variables.
func NewDSU(n int) *DSU {
	return &DSU{
		parent: intRange(n),
		parity: make([]int8, n),
		rank:   make([]int, n),
	}
}

// Find returns the component root and the x XOR root parity.
//
// The traversal is iterative so very long constraint chains cannot overflow
// the goroutine stack.
func (d *DSU) Find(x int) (root int, px int8) {
	path := []int{}
	cur := x
	for {
		path = append(path, cur)
		parent := d.parent[cur]
		if parent == cur {
			root = cur
			break
		}
		cur = parent
	}
	// Path compression from the root downwards. Before node i is rewritten,
	// d.parity[node] still holds the original edge parity, which lets the
	// running XOR accumulate the node-to-root parity in place.
	acc := int8(0)
	for i := len(path) - 1; i >= 0; i-- {
		node := path[i]
		acc = d.parity[node] ^ acc
		d.parent[node] = root
		d.parity[node] = acc
	}
	return root, acc
}

// Union records constraint x XOR y = p. The boolean reports consistency with
// previously recorded constraints.
func (d *DSU) Union(x, y int, p int8) bool {
	rx, px := d.Find(x)
	ry, py := d.Find(y)
	if rx == ry {
		return px^py == p
	}
	// Attach the smaller-rank root under the larger-rank root. When swapping
	// x and y the constraint parity p is unchanged (XOR is symmetric).
	if d.rank[rx] > d.rank[ry] {
		rx, ry = ry, rx
		px, py = py, px
	}
	d.parent[rx] = ry
	d.parity[rx] = p ^ px ^ py
	if d.rank[rx] == d.rank[ry] {
		d.rank[ry]++
	}
	return true
}

// Components assigns deterministic component ids (ordered by smallest member).
func (d *DSU) Components() []int {
	n := len(d.parent)
	ids := make([]int, n)
	next := 0
	seen := make(map[int]int, n)
	for i := 0; i < n; i++ {
		root, _ := d.Find(i)
		id, ok := seen[root]
		if !ok {
			id = next
			seen[root] = id
			next++
		}
		ids[i] = id
	}
	return ids
}

func intRange(n int) []int {
	s := make([]int, n)
	for i := range s {
		s[i] = i
	}
	return s
}

package phasing

import "testing"

func TestDSUParityConstraints(t *testing.T) {
	d := NewDSU(5)
	if !d.Union(0, 1, 1) {
		t.Fatal("first union should be accepted")
	}
	if !d.Union(1, 2, 1) {
		t.Fatal("second union should be accepted")
	}
	// 0 XOR 1 = 1 and 1 XOR 2 = 1 imply 0 XOR 2 = 0.
	if d.Union(0, 2, 1) {
		t.Fatal("expected parity contradiction for 0 XOR 2 = 1")
	}
	_, p0 := d.Find(0)
	_, p2 := d.Find(2)
	if p0^p2 != 0 {
		t.Fatalf("got parity(0)=%d parity(2)=%d, want equal", p0, p2)
	}
	ids := d.Components()
	if ids[0] != ids[1] || ids[1] != ids[2] {
		t.Fatal("0,1,2 should be in one component")
	}
	if ids[3] == ids[4] || ids[3] == ids[0] {
		t.Fatal("isolated nodes must get distinct components")
	}
}

func TestDSUPathCompression(t *testing.T) {
	d := NewDSU(64)
	for i := 0; i < 63; i++ {
		if !d.Union(i, i+1, int8(i%2)) {
			t.Fatalf("union %d failed", i)
		}
	}
	var want int8
	for i := 0; i < 63; i++ {
		want ^= int8(i % 2)
	}
	root0, p0 := d.Find(0)
	root63, p63 := d.Find(63)
	if root0 != root63 {
		t.Fatal("all nodes should share a root")
	}
	if p0^p63 != want {
		t.Fatalf("end-to-end parity = %d, want %d", p0^p63, want)
	}
	// After compression every node points directly at the root.
	for i := range d.parent {
		if d.parent[i] != root0 {
			t.Fatalf("node %d not compressed to root", i)
		}
	}
}

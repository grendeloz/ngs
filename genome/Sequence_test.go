package genome

import (
	"testing"
)

func TestSubSequence(t *testing.T) {
	s1 := NewSequence()
	s1.Name = `>chrJP | my test seq`
	s1.Sequence = `ACGTTGCA`

	e1 := s1.Sequence
	g1, err := s1.SubSequence(1, 0)
	if err != nil {
		t.Fatalf(`SubSequence 1,0 should have worked but failed: %v`, err)
	}
	if e1 != g1 {
		t.Fatalf(`SubSequence special case failed - should be %v but is %v`, e1, g1)
	}

	e2 := s1.Sequence
	g2, err := s1.SubSequence(1, 8)
	if err != nil {
		t.Fatalf(`SubSequence 1,8 should have worked but failed: %v`, err)
	}
	if e2 != g2 {
		t.Fatalf(`SubSequence 1,8 failed - should be %v but is %v`, e2, g2)
	}

	e3 := `CGTT`
	g3, err := s1.SubSequence(2, 5)
	if err != nil {
		t.Fatalf(`SubSequence 2,5 should have worked but failed: %v`, err)
	}
	if e3 != g3 {
		t.Fatalf(`SubSequence 2,5 failed - should be %v but is %v`, e3, g3)
	}

	e4 := `AC`
	g4, err := s1.SubSequence(1, 2)
	if err != nil {
		t.Fatalf(`SubSequence 1,2 should have worked but failed: %v`, err)
	}
	if e4 != g4 {
		t.Fatalf(`SubSequence 2,5 failed - should be %v but is %v`, e4, g4)
	}

	// Test the various error modes.
	g5, err := s1.SubSequence(0, 2)
	if err == nil {
		t.Fatalf(`SubSequence 0,2 should have failed but instead returned: %v`, g5)
	}
	g6, err := s1.SubSequence(2, 10)
	if err == nil {
		t.Fatalf(`SubSequence 2,10 should have failed but instead returned: %v`, g6)
	}
	g7, err := s1.SubSequence(7, 6)
	if err == nil {
		t.Fatalf(`SubSequence 7,6 should have failed but instead returned: %v`, g7)
	}
	g8, err := s1.SubSequence(9, 10)
	if err == nil {
		t.Fatalf(`SubSequence 9,10 should have failed but instead returned: %v`, g8)
	}
}

func TestWithinLimits(t *testing.T) {
	s1 := NewSequence()
	s1.Name = `>chrJP | my test seq`
	s1.Sequence = `ACGTTGCA`

	e1 := 1
	g1, ok := s1.WithinLimits(0)
	if ok != false || e1 != g1 {
		t.Fatalf(`WithinLimits(0) should have returned (%v,%v) but returned: (%v,%v)`,
			e1, false, g1, ok)
	}

	e2 := 8
	g2, ok := s1.WithinLimits(99)
	if ok != false || e2 != g2 {
		t.Fatalf(`WithinLimits(99) should have returned (%v,%v) but returned: (%v,%v)`,
			e2, false, g2, ok)
	}

	e3 := 8
	g3, ok := s1.WithinLimits(8)
	if ok != true || e3 != g3 {
		t.Fatalf(`WithinLimits(8) should have returned (%v,%v) but returned: (%v,%v)`,
			e3, false, g3, ok)
	}

	e4 := 1
	g4, ok := s1.WithinLimits(1)
	if ok != true || e4 != g4 {
		t.Fatalf(`WithinLimits(1) should have returned (%v,%v) but returned: (%v,%v)`,
			e4, false, g4, ok)
	}

	e5 := 1
	g5, ok := s1.WithinLimits(-1)
	if ok != false || e5 != g5 {
		t.Fatalf(`WithinLimits(-1) should have returned (%v,%v) but returned: (%v,%v)`,
			e5, false, g5, ok)
	}

	e6 := 7
	g6, ok := s1.WithinLimits(7)
	if ok != true || e6 != g6 {
		t.Fatalf(`WithinLimits(1) should have returned (%v,%v) but returned: (%v,%v)`,
			e6, false, g6, ok)
	}
}

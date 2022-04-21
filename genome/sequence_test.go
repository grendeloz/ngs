package genome

import (
	"testing"
)

func TestSequence(t *testing.T) {
	name := "testing"
	genome := NewGenome(name)

	file := "testdata/GRCh37_test.fa.gz"
	err := genome.AddFastaFile(file)
	if err != nil {
		t.Fatalf(`*Genome.AddFastaFile on %s failed: %v`, file, err)
	}

	e1 := 27
	g1 := len(genome.Sequences)
	if e1 != g1 {
		t.Fatalf(`Genome sequence count should be %d but is %d`, e1, g1)
	}

	seq1Name := `chr21`
	seq1, err := genome.GetSequence(seq1Name)
	if err != nil {
		t.Fatalf(`GetSequence on %s failed: %v`, seq1Name, err)
	}

	e2 := `>chr21 | 9450000 leading bases deleted (135000 lines)`
	g2 := seq1.Header
	if e2 != g2 {
		t.Fatalf(`seq Info incorrect - should be %v but is %v`, e2, g2)
	}

	e3 := `chr21`
	g3 := seq1.Name
	if e3 != g3 {
		t.Fatalf(`seq Name incorrect - should be %v but is %v`, e3, g3)
	}

	e4 := `9450000 leading bases deleted (135000 lines)`
	g4 := seq1.Info
	if e4 != g4 {
		t.Fatalf(`seq Info incorrect - should be %v but is %v`, e4, g4)
	}
}

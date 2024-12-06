package genome

import (
	"testing"
)

func TestGenomeAddFastaFile(t *testing.T) {
	name := "testing"
	genome := NewGenome(name)

	// Check defaulting
	if genome.Name != name {
		t.Fatalf(`Genome name should be %s but is %s`, name, genome.Name)
	}

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

	e2 := `d900cf9254cc50cbce326163f78acebc`
	g2, err := genome.Sequences[0].FastaFile.MD5()
	if err != nil {
		t.Fatalf(`*MD5() onn %s failed: %v`, file, err)
	}
	if e2 != g2 {
		t.Fatalf(`seq 0 FastaFile.MD5 incorrect - should be %v but is %v`, e2, g2)
	}

	e3 := `testdata/GRCh37_test.fa.gz`
	g3 := genome.Sequences[0].FastaFile.Filepath
	if e3 != g3 {
		t.Fatalf(`seq 0 FastaFile.Filepath incorrect - should be %v but is %v`, e3, g3)
	}
}

package genome

import (
	"testing"
)

var fa1 = `>chr1
ACGTCCAGCC
GACTCGg
AGCGACGA
>chr2
ACGTCCAGCCGACTCgGCGACGA
>chr3
CGTCCAGCC
GACTCGg
AGCGACGA
`

// TestDefaults tests that a new Reader has expected defaults
func TestFastaDefaults(t *testing.T) {
	name := "testing"
	genome := NewGenome(name)

	// Check defaulting
	if genome.Name != name {
		t.Fatalf(`Genome name should be %s but is %s`, name, genome.Name)
	}

	e1 := 0
	g1 := len(genome.Sequences)
	if e1 != g1 {
		t.Fatalf(`Genome sequence count should be %d but is %d`, e1, g1)
	}
}

// TestOpenFastaFile
func TestOpenFastaFile(t *testing.T) {
	file := "testdata/GRCh37_test.fa.gz"

	// Check
	seqs, err := OpenFastaFile(file)
	if err != nil {
		t.Fatalf(`OpenFastaFile on %s failed: %v`, file, err)
	}

	e1 := 27
	g1 := len(seqs)
	if e1 != g1 {
		t.Fatalf(`Genome sequence count should be %d but is %d`, e1, g1)
	}

	e2 := 16569
	g2 := seqs[26].Length()
	if e2 != g2 {
		t.Fatalf(`seq 26 (chrMT) length incorrect - should be %d but is %d`, e2, g2)
	}

	// Test that incorrect value doesn't succeed
	e3 := 16568
	g3 := seqs[26].Length()
	if e3 == g3 {
		t.Fatalf(`seq 26 (chrMT) length incorrect - should not be %d but is %d`, e3, g3)
	}

	e7 := `>chr1 | 70000 leading bases deleted`
	g7 := seqs[0].Header
	if e7 != g7 {
		t.Fatalf(`seq 0 Header incorrect - should be %v but is %v`, e7, g7)
	}

	e4 := 700
	g4 := seqs[0].Length()
	if e4 != g4 {
		t.Fatalf(`seq 0 (chr1) length incorrect - should be %d but is %d`, e4, g4)
	}

	e5 := `AGTTTTAGAT`
	g5 := seqs[0].Sequence[0:10]
	if e5 != g5 {
		t.Fatalf(`seq 0 sequence incorrect - should be %v but is %v`, e5, g5)
	}

	e6 := `TCCATCTATG`
	g6 := seqs[0].Sequence[690:700]
	if e6 != g6 {
		t.Fatalf(`seq 0 sequence incorrect - should be %v but is %v`, e6, g6)
	}

	e11 := `>chr21 | 9450000 leading bases deleted (135000 lines)`
	g11 := seqs[20].Header
	if e11 != g11 {
		t.Fatalf(`seq 20 Header incorrect - should be %v but is %v`, e11, g11)
	}

	e10 := `CCAGGTTCAA`
	g10 := seqs[20].Sequence[700:710]
	if e10 != g10 {
		t.Fatalf(`seq 20 sequence incorrect - should be %v but is %v`, e10, g10)
	}

	e20 := `>chrMT`
	g20 := seqs[26].Header
	if e20 != g20 {
		t.Fatalf(`seq 26 Header incorrect - should be %v but is %v`, e20, g20)
	}

	e30 := `d900cf9254cc50cbce326163f78acebc`
	g30 := seqs[26].FastaFile.MD5
	if e30 != g30 {
		t.Fatalf(`seq 26 FastaFile.MD5 incorrect - should be %v but is %v`, e30, g30)
	}

	e31 := `testdata/GRCh37_test.fa.gz`
	g31 := seqs[26].FastaFile.Filepath
	if e31 != g31 {
		t.Fatalf(`seq 26 FastaFile.Filepath incorrect - should be %v but is %v`, e31, g31)
	}

	e32 := `d900cf9254cc50cbce326163f78acebc`
	g32 := seqs[0].FastaFile.MD5
	if e32 != g32 {
		t.Fatalf(`seq 0 FastaFile.MD5 incorrect - should be %v but is %v`, e32, g32)
	}

	e33 := `testdata/GRCh37_test.fa.gz`
	g33 := seqs[0].FastaFile.Filepath
	if e33 != g33 {
		t.Fatalf(`seq 0 FastaFile.Filepath incorrect - should be %v but is %v`, e33, g33)
	}
}

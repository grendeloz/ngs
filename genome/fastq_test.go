package genome

import (
	"testing"
)

// TestParseFastqFile
func TestParseFastqFile(t *testing.T) {
	type test struct {
		Id   string
		Seq  string
		Qual string
	}

	tests := []test{
		{"@read1",
			"ACGTCCAGCCACGTCCAGCCGACTCGGCGA",
			"ABCDEFGHIJLKMNOPQRSTUVWXYZ1234"},
		{"@read2",
			"CGTCCAGCCACGTCCAGCCGACTCGGCGAA",
			"BCDEFGHIJLKMNOPQRSTUVWXYZ12345"},
		{"@read3",
			"GTCCAGCCACGTCCAGCCGACTCGGCGAAC",
			"CDEFGHIJLKMNOPQRSTUVWXYZ123456"},
	}

	file := "testdata/test1.fq"

	// Check that FastqFile initialises
	ff, err := NewFastqFile(file)
	if err != nil {
		t.Fatalf(`NewFastqFile on %s failed: %v`, file, err)
	}

	// Check header
	e1 := 1
	g1 := len(ff.Headers)
	if e1 != g1 {
		t.Fatalf(`header line count incorrect - should be %d but is %d`, e1, g1)
	}
	e2 := "# header line"
	g2 := ff.Headers[0]
	if e2 != g2 {
		t.Fatalf(`header line incorrect - should be [%s] but is [%s]`, e2, g2)
	}

	// Check read 1
	for i, tst := range tests {
		found := ff.Next()
		if !found {
			t.Fatalf(`read %d not found by Next()`, i)
		}
		rec := ff.Record()
		if tst.Id != rec.Id {
			t.Fatalf(`read %d Id incorrect - expected %s got %s`,
				i, tst.Id, rec.Id)
		}
		if tst.Seq != rec.Seq {
			t.Fatalf(`read %d Seq incorrect - expected %s got %s`,
				i, tst.Seq, rec.Seq)
		}
		if tst.Qual != rec.Qual {
			t.Fatalf(`read %d Qual incorrect - expected %s got %s`,
				i, tst.Qual, rec.Qual)
		}
	}
}

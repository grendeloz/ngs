package genome

import (
	"testing"
)

func TestOpenFastqFile(t *testing.T) {
	type test struct {
		Id   string
		Seq  string
		Qual string
	}

	// Records that are found in testdata/test1.fq
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
	ff, err := OpenFastqFile(file)
	if err != nil {
		t.Fatalf(`OpenFastqFile on %s failed: %v`, file, err)
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
		rec, err := ff.Next()
		if err != nil {
			t.Fatalf(`Next() threw an unexpected error: %v`, err)
		}
		if rec == nil {
			t.Fatalf(`test read %d not found by Next()`, i)
		}

		if tst.Id != rec.Id {
			t.Fatalf(`read %d Id incorrect - expected %s got %s`,
				i, tst.Id, rec.Id)
		}
		seq := string(rec.Bases)
		if tst.Seq != seq {
			t.Fatalf(`read %d Seq incorrect - expected %s got %s`,
				i, tst.Seq, seq)
		}
		qual := string(rec.Qualities)
		if tst.Qual != qual {
			t.Fatalf(`read %d Qual incorrect - expected %s got %s`,
				i, tst.Qual, qual)
		}
	}
}

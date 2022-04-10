// package header processes comments from the start of a text file.
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

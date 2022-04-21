package genome

import (
	"log"
	"regexp"
	"strings"
)

type Sequence struct {
	// Header is the complete string from the '>' line that is at the
	// start of each sequence in FASTA. The Header is often a composite
	// |-separated list of elements and so may not be directly useful.
	Header string

	// Regardless of the format of the Header line, the sequence Name is
	// always the first string.
	Name string

	// If the Header line contains more that just the sequence Name,
	// this string holds the rest of the Header line. Any spaces or |
	// between the Name and the Info have been removed as have any
	// trailing spaces.
	Info string

	// The base sequence.
	Sequence string

	// A link to the FASTA file object from which this Sequence was
	// generated.
	FastaFile *FastaFile
}

// NewSequence takes the header line from the FASTA and returns a new
// Sequence object.
func NewSequence(header string) *Sequence {
	seq := &Sequence{Header: header}
	// Parse the header
	header = strings.TrimLeft(header, " >")
	header = strings.TrimSpace(header)

	re := regexp.MustCompile(`([^ \t]+)([ |]*)(.*)`)
	matches := re.FindStringSubmatch(header)
	switch len(matches) {
	case 2:
		seq.Name = matches[1]
	case 4:
		seq.Name = matches[1]
		seq.Info = matches[3]
	default:
		log.Fatalf("Header failed parsing: %s", header)
	}

	return seq
}

func (s *Sequence) Length() int {
	return len(s.Sequence)
}

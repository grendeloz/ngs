package genome

import (
	"regexp"
	"strings"
)

var faParserRex *regexp.Regexp = regexp.MustCompile(`([^ \t]+)([ |]*)(.*)`)

type FastaRec struct {
	// Header is the complete string from the '>' line that is at the
	// start of each sequence in FASTA. The Header is often a composite
	// |-separated list of elements and so may not be directly useful.
	Header string

	// Regardless of the format of the Header line, the sequence Name is
	// always the start of the Header string up to the first space or
	// pipe symbol "|".
	Name string

	// If the Header line contains more that just the sequence Name,
	// this string holds the rest of the Header line. Any spaces or |
	// between the Name and the Info have been removed as have any
	// trailing spaces.
	Info string

	// Base sequence
	Sequence string

	// FASTA file object that contained this record
	FastaFile *FastaFile
}

// NewFastaRec takes the header line from the FASTA and returns a new
// object.
func NewFastaRec(header string) *FastaRec {
	r := &FastaRec{Header: header}
	// Parse the header
	header = strings.TrimLeft(header, " >")
	header = strings.TrimSpace(header)

	matches := faParserRex.FindStringSubmatch(header)
	switch len(matches) {
	case 2:
		r.Name = matches[1]
	case 4:
		r.Name = matches[1]
		r.Info = matches[3]
	default:
		// weird Header so do nothing
	}

	return r
}

// Length is the length of the sequence.
func (r *FastaRec) Length() int {
	return len(r.Sequence)
}

package genome

import (
	"fmt"
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

// Length is the length of the sequence.
func (s *Sequence) Length() int {
	return len(s.Sequence)
}

// WithinLimits checks whether the supplied integer is within the
// start-end limits of a sequence. If it is within the limits, the input
// integer is returned and the bool is set to true. If the input is
// before the beginning or beyond the end of the sequence, the relevant
// limit is returned and the bool is set to false. It might
// be used like this:
//
//   t := 16569
//   l, ok := s.WithinLimit(t)
//   if ! ok {
//       fmt.Println("Position %d is beyond sequence limits: %d", t, l)
//   }
//
// WithinLimits assumes sequences are 1-based closed intervals, i.e.
// a sequence with 10 bases starts at 1 and ends at 10. For this 10bp
// sequence, WithinLimits(10) is true and WithinLimits(0) and
// WithinLimits(11) are both false.
func (s *Sequence) WithinLimits(t int) (int, bool) {
	switch {
	case t < 1:
		return 1, false
	case t > len(s.Sequence):
		return len(s.Sequence), false
	default:
		return t, true
	}
}

// SubSequence returns a subsequence of Sequence.Sequence. The start and
// end values form a 1-based closed interval, i.e. the start and end
// bases are both included in the subsequence. Setting end to 0 is a
// special case - it signals that you want all sequence from the start
// position to the end of the sequence.
//
// Given a sequence ACGTTGCA, this list shows what would be returned if
// we were to call SubSequence with the start,end values:
//
//   1,4 : ACGT - 4-base sequence from the 1st to 4th base inclusive.
//   0,4 : error, start cannot be less than 1
//   1,8 : ACGTTGCA - returns the entire sequence
//   1,0 : ACGTTGCA - returns the entire sequence, special case
//   1,9 : error - end cannot be > length(sequence)
//   4,3 : error - start cannot be > end
//   5,5 : T - returns a single base
func (s *Sequence) SubSequence(start, end int) (string, error) {
	// Note that golang substrings are 0-based half-open intervals, i.e.
	// the start is included but the end is the first base past the end of
	// the subsequence. So 1,4 returns the 2nd and 3rd bases. It is on
	// us to do the conversion.
	switch {
	case start < 1:
		return "", fmt.Errorf("SubSequence: start cannot be less than 1: %d", start)
	case end > len(s.Sequence):
		return "", fmt.Errorf("SubSequence: end cannot be beyond the end of the sequence: %d", end)
	case start > len(s.Sequence):
		return "", fmt.Errorf("SubSequence: start cannot be beyond the end of the sequence: %d", start)
	case end == 0:
		// must come before start>end case
		end = len(s.Sequence)
	case start > end:
		return "", fmt.Errorf("SubSequence: start cannot be > end: %d", start)
	}
	// Convert from 1-based closed sequence coords
	//           to 0-based half-open go string coords ...
	return s.Sequence[start-1 : end], nil
}

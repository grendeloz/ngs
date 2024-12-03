package genome

import (
	"fmt"
	"strings"
)

type FastqRec struct {
	// Id is the string that appears as the first line of a 4-line record.
	// In a FASTQ file, the Id must start with a "@" character but it is
	// removed here. The Id may consist of a single "word" or it may
	// contain spaces which separate multiple words. The third line of a
	// record must start with a "+" character and then may be either
	// empty or repeat the Id (without the leading "@").
	Id string

	// The base sequence.
	Bases []byte

	// The base quality scores encoded as ASCII characters. It is an
	// error for the lengths of the Bases and Qualities arrays to be of
	// different lengths.
	Qualities []byte
}

// NewFastqRec returns an empty FastqRec.
func NewFastqRec() *FastqRec {
	r := FastqRec{}
	r.Bases = make([]byte, 0)
	r.Qualities = make([]byte, 0)
	return &r
}

// FastqRecFromString assumes that a string with 4 newline-terminated
// lines will be supplied.
func FastqRecFromString(s string) (*FastqRec, error) {
	r := NewFastqRec()
	lines := strings.Split(s, "\n")
	if len(lines) < 4 {
		return r, fmt.Errorf("a FASTQ record has 4 lines, input string has: %d", len(lines))
	}

	// Parse Id
	var id []byte = []byte(lines[0])
	if id[0] != '@' {
		return r, fmt.Errorf("first char of id line must be @ but is [%c]", id[0])
	}
	r.Id = string(id[1:])

	var id2 []byte = []byte(lines[2])
	if id2[0] != '+' {
		return r, fmt.Errorf("first char of second id line must be + but is [%c]", id2[0])
	}

	r.SetBasesFromString(lines[1])
	r.SetQualitiesFromString(lines[3])
	if err := r.CheckValid(); err != nil {
		return r, fmt.Errorf("read is not valid: %v", err)
	}

	return r, nil
}

func (r *FastqRec) SetBasesFromString(s string) {
	r.Bases = []byte(s)
}

func (r *FastqRec) SetQualitiesFromString(s string) {
	r.Qualities = []byte(s)
}

// CheckValid checks that a Record has an Id and that the count of Bases
// and Qualities is the same. Note that a Record with no Bases and no
// Qualities is considered valid.
func (r *FastqRec) CheckValid() error {
	if len(r.Bases) != len(r.Qualities) {
		return fmt.Errorf("base and quality score counts do not match for read: %s", r.Id)
	}
	return nil
}

// String returns a 4-line string representation of the Record with "\n"
// as the line-ending and "+" by itself for line 3.
func (r *FastqRec) String() string {
	var builder strings.Builder
	var pieces []string

	pieces = append(pieces, "@", r.Id, "\n",
		string(r.Bases), "\n",
		"+\n",
		string(r.Qualities), "\n")

	// We are ignoring errors thrown by WriteString(). Naughty.
	for _, p := range pieces {
		builder.WriteString(p)
	}
	return builder.String()

}

package genome

import (
	"fmt"
	"strings"
)

type Read struct {
	// Id is the string that appears as the first line of a 4-line read.
	// In a file, the Id must start with a "@" character which is
	// removed for the Id field here in Read. The Id may consist of a
	// single "word" or it may contain spaces which separate multiple
	// words. The third line of a read must start with a "+" character
	// and then may either be empty or repeat the Id (without the leading
	// "@").
	Id string

	// The base sequence.
	Bases []byte

	// The base quality scores encoded as ASCII characters. It is an
	// error for the lengths of the Bases and Qualities arrays to be of
	// different lengths.
	Qualities []byte
}

// NewRead returns an empty Read.
func NewRead() *Read {
	r := Read{}
	r.Bases = make([]byte, 0)
	r.Qualities = make([]byte, 0)
	return &r
}

// NewReadFromString assumes that a string with 4 newline-terminated
// lines will be supplied.
func NewReadFromString(s string) (*Read, error) {
	r := &Read{}
	lines := strings.Split(s, "\n")
	if len(lines) < 4 {
		return r, fmt.Errorf("a read has 4 lines, input string has: %d", len(lines))
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

func (r *Read) SetBasesFromString(s string) {
	r.Bases = []byte(s)
}

func (r *Read) SetQualitiesFromString(s string) {
	r.Qualities = []byte(s)
}

// CheckValid checks that a Read has an Id and that the count of Bases
// and Qualities is the same. Note that a Read with no Bases and no
// Qualities is considered valid.
func (r *Read) CheckValid() error {
	if len(r.Bases) != len(r.Qualities) {
		return fmt.Errorf("base and quality score counts do not match for read: %s", r.Id)
	}
	return nil
}

// String returns a 4-line string representation of the Read with "\n"
// as the line-ending and "+" by itself for line 3.
func (r *Read) String() string {
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

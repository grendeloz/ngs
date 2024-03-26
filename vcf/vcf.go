// Package vcf is a lightweight reader/writer for genomics VCF files.
// It is based on the VCF 4.3 spec foudn at:
// https://samtools.github.io/hts-specs/VCFv4.3.pdf
package vcf

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"os"
	"regexp"
	"strings"
)

// Errors
var (
	ErrNoVcfMeta = errors.New("No meta lines - not compliant with VCF spec")
)

// Patterns for text parsing
var metaRx = regexp.MustCompile(`^##`)
var headRx = regexp.MustCompile(`^#`)

// A Selector defines a selection operation, the type of thing to be
// selected, and the pattern (regex) to use for selection.
type Vcf struct {
	Meta    Meta
	Header  Header
	Records Records
}

func (v Vcf) String() string {
	return v.Meta.String() +
		v.Header.String() +
		v.Records.String()
}

type Meta struct {
	OrigStr string // string as read from file
}

type Header struct {
	OrigStr string // string as read from file
}

type Records struct {
	OrigStr string // string as read from file
}

func NewVcf() *Vcf {
	return &Vcf{Meta: NewMeta(),
		Header:  NewHeader(),
		Records: NewRecords()}
}

// Clone creates a deep copy of a Vcf, i.e. the new Vcf shares no
// pointers with the original Vcf. After calling Clone you can change
// the original Vcf or the copy without any concern that the change
// will affect the other Vcf. Obviously your memory use doubles. This is
// useful when you want to filter or change a VCF file because it lets
// you copy the original and change the records while leaving the header
// and column names intact.
func (v *Vcf) Clone() *Vcf {
	nv := NewVcf()
	nv.Meta = v.Meta.Clone()
	nv.Header = v.Header.Clone()
	nv.Records = v.Records.Clone()
	return nv
}

// NewFromFile reads from a file and returns a pointer to a Vcf.
func NewFromFile(file string) (*Vcf, error) {
	// Open file
	ff, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer ff.Close()

	// We need to define this before we handle gzip
	var scanner *bufio.Scanner

	// Based on file extension, handle gzip files
	found, err := regexp.MatchString(`\.[gG][zZ]$`, file)
	if err != nil {
		return nil, fmt.Errorf("NewFromFile: error matching gzip file pattern against %s: %w", file, err)
	}
	if found {
		// For gzip files, put a gzip.Reader into the chain
		reader, err := gzip.NewReader(ff)
		if err != nil {
			return nil, fmt.Errorf("NewFromFile: error opening gzip file %s: %w", file, err)
		}
		defer reader.Close()
		scanner = bufio.NewScanner(reader)
	} else {
		// For non gzip files, go straight to bufio.Reader
		scanner = bufio.NewScanner(ff)
	}

	vcf, err := newFromScanner(scanner)
	if err != nil {
		return vcf, fmt.Errorf("NewFromFile: error scanning: %w", err)
	}
	return vcf, nil
}

// newFromScanner reads from a *bufio.Scanner and returns a pointer
// to a Vcf. It is an alternative to NewFromFile and is useful when
// you have a VCF file as a block of text in memory.
func newFromScanner(scanner *bufio.Scanner) (*Vcf, error) {
	vcf := NewVcf()

	// Unnecessary but explicit
	scanner.Split(bufio.ScanLines)

	// Let's do string concatenation the fast way
	var mb, hb, rb strings.Builder

	// Read the file
	for scanner.Scan() {
		//line := strings.TrimSuffix(scanner.Text(), "\n")
		if metaRx.MatchString(line) {
			mb.WriteString(line)
		} else if headRx.MatchString(line) {
			hb.WriteString(line)
		} else {
			rb.WriteString(line)
		}
	}

	// If there are no Meta lines then it can't be a VCF because the
	// fileformat= meta line as the first line is mandatory.
	if len(mb.String()) == 0 {
		return nil, ErrNoVcfMeta
	}

	// If the very first line is not the gff3 identifier then we exit
	// immediately. Unfortunately, some folks (Ensembl at least) seem to
	// put an arbitrary number of spaces in as separator so we will need
	// to use a pattern rather than a simple string equality test.
	ok, err := regexp.Match(`^##fileformat=VCFv`, []byte(mb.String()))
	if err != nil {
		return nil, fmt.Errorf("newFromScanner: error matching fileformat line: %w", err)
	}
	if !ok {
		return nil, fmt.Errorf("newFromScanner: mandatory fileformat= line missing", gff3.Header[0])
	}

	vcf.Meta.OrigStr = mb.String()
	vcf.Header.OrigStr = hb.String()
	vcf.Records.OrigStr = rb.String()

	return vcf, nil
}

func (v *Vcf) Write(file string) error {
	f, err := os.Create(file)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	// TO DO
	// this all needs to change because this just writes out the
	// original string which is obviously not what we want.

	_, err = w.WriteString(v.Meta.OrigStr + v.Header.OrigStr + v.Records.OrigStr)
	if err != nil {
		return err
	}

	return nil
}

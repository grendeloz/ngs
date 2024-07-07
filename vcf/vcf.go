// Package vcf is a lightweight reader/writer for genomics VCF files.
// It is based on the VCF 4.3 specification found at:
// https://samtools.github.io/hts-specs/VCFv4.3.pdf. Note that vcf
// does not attempt to enforce the check for the validity of a VCF file
// and the lines within. vcf is mosty concerned with parsing the
// contents of the file and anything beyond that is down to the user.
// While we consulted the VCF 4.3 specification while creating vcf, it
// will probably successfully parse older and newer VCF files.
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
	ErrNoVcfMeta    = errors.New("No meta lines - not compliant with VCF spec")
	ErrNoFileformat = errors.New("Does not start with valid fileformat Meta line")
)

// Patterns for text parsing.
var fileformatRx = regexp.MustCompile(`^##fileformat=version ([^\s]+)$`)
var metaStructuredRx = regexp.MustCompile(`^##([^=]+)=<(.+)>`)
var metaUnstructuredRx = regexp.MustCompile(`^##([^=]+)=([^<].*)`)
var headRx = regexp.MustCompile(`^#[^#]`)
var gzipRx = regexp.MustCompile(`\.[gG][zZ]$`)

func IsGzip(filename string) bool {
	return gzipRx.MatchString(filename)
}

func IsFileformatMeta(line string) bool {
	return fileformatRx.MatchString(line)
}

func IsMetaUnstructured(line string) bool {
	return fileformatRx.MatchString(line)
}

// Vcf holds information parsed from a VCF file. If ReadAll() is used,
// then the VCF Records are in Records. If Next() is used then Records
// starts empty and stays empty unless you manually add Records to it.
type Vcf struct {
	Meta       *Meta
	Header     *Header
	Records    []*Record
	Fileformat string
	mOrigStr   string
	hOrigStr   string
	rOrigStr   string
}

// String does what you would expect. It is a simple way to get a text
// representation of the whole VCF so you can write it out but it does
// take up a considerable amount of memory to hold the text
// representation.
func (v *Vcf) String() string {
	s := v.Meta.String() + "\n" +
		v.Header.String() + "\n"
	for _, r := range v.Records {
		s = s + r.String() + "\n"
	}
	return s
}

func NewVcf() *Vcf {
	return &Vcf{Meta: NewMeta(),
		Header:  NewHeader(),
		Records: make([]*Record, 0, 100)}
}

/*
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
*/

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

	// Based on file extension, handle gzip files. For gzip files,
	// put a gzip.Reader into the chain. For non-gzip files, go
	// straight to a bufio.Reader
	if IsGzip(file) {
		reader, err := gzip.NewReader(ff)
		if err != nil {
			return nil, fmt.Errorf("NewFromFile: error opening gzip file %s: %w", file, err)
		}
		defer reader.Close()
		scanner = bufio.NewScanner(reader)
	} else {
		scanner = bufio.NewScanner(ff)
	}

	// Parse Meta and Header lines.
	vcf, err := newFromScanner(scanner)
	if err != nil {
		return vcf, fmt.Errorf("NewFromFile: error scanning: %w", err)
	}
	return vcf, nil
}

// newFromScanner reads from a *bufio.Scanner and returns a pointer
// to a Vcf. Because it reads from a Scanner, it work equally well
// with Scanners against files or Scanners against strings in memory.
// It is used within NewFromFile().
func newFromScanner(scanner *bufio.Scanner) (*Vcf, error) {
	vcf := NewVcf()

	// Unnecessary but explicit
	scanner.Split(bufio.ScanLines)

	// Let's do string concatenation the fast way
	var mb, hb, rb strings.Builder

	// Read everything except the records. Structure must be:
	// - a fileformat Meta line
	// - zero or more structured or unstructured Meta lines
	// - a Header line
	var line string

	scanner.Scan()
	line = scanner.Text()
	if IsMetaUnstructured(line) && IsFileformatMeta(line) {
		m := fileformatRx.FindStringSubmatch(line)
		vcf.Fileformat = m[1]
	} else {
		return nil, ErrNoFileformat
	}

	for scanner.Scan() {
		//line := strings.TrimSuffix(scanner.Text(), "\n")
		line = scanner.Text()
		if metaStructuredRx.MatchString(line) {
			mb.WriteString(line)
		} else if metaUnstructuredRx.MatchString(line) {
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
		return nil, fmt.Errorf("newFromScanner: mandatory fileformat= first line missing")
	}

	vcf.mOrigStr = mb.String()
	vcf.hOrigStr = hb.String()
	vcf.rOrigStr = rb.String()

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

	_, err = w.WriteString(v.mOrigStr + v.hOrigStr + v.rOrigStr)
	if err != nil {
		return err
	}

	return nil
}

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
var fileformatRx = regexp.MustCompile(`^##fileformat=(.*)$`)
var metaStructuredRx = regexp.MustCompile(`^##([^=]+)=<(.+)>$`)
var metaUnstructuredRx = regexp.MustCompile(`^##([^=]+)=([^<].*)$`)
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
	Samples    []string
	Fileformat string
	mOrigStr   string
	hOrigStr   string
	rOrigStr   string
}

// String does what you would expect. It is a simple way to get a text
// representation of the whole VCF so you can write it out but be aware
// that it may take up a considerable amount of memory to create the
// string representation while also holding the data structure.
func (v *Vcf) String() string {
	// lazy lazy lazy
	s := v.mOrigStr + v.hOrigStr

	// TO DO
	// Make this work and then use the same logic in Write() so we
	// stream the string representation to disk rather than creating a
	// full string representation in memory and then writing it out.
	//s := v.Meta.String() + v.Header.String()
	//for _, r := range v.Records {
	//	s = s + r.String() + "\n"
	//}

	return s
}

func NewVcf() *Vcf {
	return &Vcf{Meta: NewMeta(),
		Header:  NewHeader(),
		Records: make([]*Record, 0, 100),
		Samples: make([]string, 0)}
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
	//fmt.Printf("first line read: %v\n", line)
	if IsFileformatMeta(line) {
		m := fileformatRx.FindStringSubmatch(line)
		vcf.Fileformat = m[1]

		mb.WriteString(line)
		mb.WriteByte('\n')
	} else {
		return nil, ErrNoFileformat
	}

	var mUn, mSt, h, r int
	mUn = 1 // fileformat line

	for scanner.Scan() {
		//line := strings.TrimSuffix(scanner.Text(), "\n")
		line = scanner.Text()
		if metaStructuredRx.MatchString(line) {
			mSt++
			mb.WriteString(line)
			mb.WriteByte('\n')
		} else if metaUnstructuredRx.MatchString(line) {
			mUn++
			mb.WriteString(line)
			mb.WriteByte('\n')
		} else if headRx.MatchString(line) {
			h++
			hb.WriteString(line)
			hb.WriteByte('\n')
		} else {
			r++
			rb.WriteString(line)
			rb.WriteByte('\n')
			r, err := RecordFromString(line)
			if err != nil {
				return vcf, fmt.Errorf("problem parsing record %s: %w", line, err)
			}
			vcf.Records = append(vcf.Records, r)
		}
	}
	fmt.Printf("line counts: mUn:%d mSt:%d Header:%d Records:%d\n", mUn, mSt, h, r)

	// If there are no Meta lines then it can't be a VCF because the
	// fileformat= meta line as the first line is mandatory.
	if len(mb.String()) == 0 {
		return nil, ErrNoVcfMeta
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

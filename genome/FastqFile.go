package genome

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"regexp"
	"strings"
)

// Pattern for header lines
var fqHeaderRex *regexp.Regexp = regexp.MustCompile(`^#(.*)$`)

// FastqFile
type FastqFile struct {
	Filepath  string
	Headers   []string
	scanner   *bufio.Scanner // used in Next()
	recCtr    int
	md5       string
	nextRecId string
	EOF       bool
}

// OpenFastqFile opens a FASTQ file and prepares it for reading.
// It will handle gzipped files as long as they have a .gz extension.
func OpenFastqFile(file string) (*FastqFile, error) {
	// As a side effect of reading the FASTQ
	fastq := &FastqFile{Filepath: file}
	fastq.Headers = make([]string, 0)

	// Do NOT defer close of readers - we want them to stay open
	// and be passed around in FastqFile.
	f, err := os.Open(file)
	if err != nil {
		return fastq, err
	}

	// We need to define this before we handle gzip
	//var scanner *bufio.Scanner

	// Based on file extension, handle gzip files
	found, err := regexp.MatchString(`\.[gG][zZ]$`, file)
	if err != nil {
		return fastq, fmt.Errorf("error matching gzip file pattern: %v", err)
	}
	if found {
		// For gzip files, put a gzip.Reader into the chain
		reader, err := gzip.NewReader(f)
		if err != nil {
			return fastq, fmt.Errorf("unable to open gzip file %v: %w", file, err)
		}
		fastq.scanner = bufio.NewScanner(reader)
	} else {
		// For non gzip files, go straight to bufio.Reader
		fastq.scanner = bufio.NewScanner(f)
	}

	// Unnecessary but explicit
	fastq.scanner.Split(bufio.ScanLines)

	// Read header
	for fastq.scanner.Scan() {
		line := strings.TrimSuffix(fastq.scanner.Text(), "\n")
		if fqHeaderRex.MatchString(line) {
			fastq.Headers = append(fastq.Headers, line)
		} else {
			// Found Id line of the first record so save and return
			fastq.nextRecId = line
			return fastq, nil
		}
	}

	return fastq, nil
}

// Next returns the next record from the FASTQ file. If there are no
// more records, it returns nil.
func (f *FastqFile) Next() (*FastqRec, error) {
	if f.EOF {
		return nil, nil
	}

	thisRec := NewFastqRec()
	f.recCtr++

	// First record special case - we already read the first line
	if f.nextRecId != "" {
		thisRec.Id = f.nextRecId
		f.nextRecId = ""
	} else {
		f.scanner.Scan()
		thisRec.Id = f.scanner.Text()
	}

	// Read the next 3 lines
	f.scanner.Scan()
	thisRec.Bases = []byte(f.scanner.Text())
	f.scanner.Scan()
	f.scanner.Scan()
	thisRec.Qualities = []byte(f.scanner.Text())

	return thisRec, nil
}

// RecordCount returns the number of records returned with Next().
func (f *FastqFile) RecordCount() int {
	return f.recCtr
}

// MD5 will return the MD5 string for the file and will calculate it on
// the first call, which is therefore slow. Subsequent calls return the
// already-calculated value.
func (f *FastqFile) MD5() (string, error) {
	if f.md5 != "" {
		return f.md5, nil
	}
	md5, err := Md5sum(f.Filepath)
	if err != nil {
		return md5, fmt.Errorf("error generating MD5 for %v: %w", f.Filepath, err)
	}
	f.md5 = md5
	return f.md5, nil
}

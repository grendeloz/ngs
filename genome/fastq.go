package genome

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"regexp"
)

// Pattern for header lines
var hrex *regexp.Regexp = regexp.MustCompile(`^#(.*)$`)

// FastqFile must be uniquely identifiable so we can check that a given
// FASTQ file matches the one used to create Genomes etc. Files can
// change location so Filepath is not enough - MD5 is not a perfect
// solution but it's a pretty handy siganture.
type FastqFile struct {
	Filepath string
	Headers  []string
	scanner  *bufio.Scanner
	rec      *FastqRec
	recCtr   int
}

type FastqRec struct {
	Id   string
	Seq  string
	Qual string
}

func NewFastqFile(file string) (*FastqFile, error) {
	fastq := &FastqFile{Filepath: file}
	fastq.Headers = make([]string, 0)
	fastq.rec = &FastqRec{}

	// Do not defer close of readers because we want them to stay open
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
		return fastq, fmt.Errorf("error matching pattern for gzip file: %v", err)
	}
	if found {
		// For gzip files, put a gzip.Reader into the chain
		reader, err := gzip.NewReader(f)
		if err != nil {
			return fastq, fmt.Errorf("unable to open gzip file %v: %v", file, err)
		}
		fastq.scanner = bufio.NewScanner(reader)
	} else {
		// For non gzip files, go straight to bufio.Reader
		fastq.scanner = bufio.NewScanner(f)
	}
	//fastq.scanner = scanner

	fastq.scanner.Split(bufio.ScanLines) // unnecessary but explicit

	// Read header
	for fastq.scanner.Scan() {
		if hrex.MatchString(fastq.scanner.Text()) {
			fastq.Headers = append(fastq.Headers, fastq.scanner.Text())
		} else {
			// we have hit the first line of the first record
			// so read the other 3 lines of the record
			fastq.rec.Id = fastq.scanner.Text()
			fastq.scanner.Scan()
			fastq.rec.Seq = fastq.scanner.Text()
			fastq.scanner.Scan()
			fastq.scanner.Scan()
			fastq.rec.Qual = fastq.scanner.Text()
			break
		}
	}

	return fastq, nil
}

// Next does not return a record but it sets up the next record for
// retrieval by Record(). It is designed to be used
func (f *FastqFile) Next() bool {
	// On first invocation, do nothing because we already have the first
	// record in place because it was read during NewFastqReader
	if f.recCtr == 0 {
		f.recCtr++
		return true
	}

	// Read First line to see if there is a next record
	found := f.scanner.Scan()
	if !found {
		return false
	}

	// Read the next 4 lines into ff.rec
	f.rec.Id = f.scanner.Text()
	f.scanner.Scan()
	f.rec.Seq = f.scanner.Text()
	f.scanner.Scan()
	f.scanner.Scan()
	f.rec.Qual = f.scanner.Text()

	return true
}

func (f *FastqFile) Record() *FastqRec {
	// Return a copy of the current record
	return &FastqRec{Id: f.rec.Id,
		Seq:  f.rec.Seq,
		Qual: f.rec.Qual}
}

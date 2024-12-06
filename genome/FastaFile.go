package genome

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"regexp"
	"strings"
)

// Pattern for header (comment) and Id lines
var faHeaderRex *regexp.Regexp = regexp.MustCompile(`^;(.*)$`)
var faIdRex *regexp.Regexp = regexp.MustCompile(`^>(.*)$`)

// FastaFile
type FastaFile struct {
	Filepath  string
	Headers   []string
	scanner   *bufio.Scanner // used in Next()
	recCtr    int
	md5       string
	nextRecId string
	EOF       bool
}

// OpenFastaFile opens a FASTA file and prepares it for reading.
// It will handle gzipped files as long as they have a .gz extension.
func OpenFastaFile(file string) (*FastaFile, error) {
	fasta := &FastaFile{Filepath: file}
	fasta.Headers = make([]string, 0)

	// Do NOT close or defer close readers - we want them to stay open
	// and be passed around in FastaFile.
	ff, err := os.Open(file)
	if err != nil {
		return fasta, err
	}

	// Based on file extension, handle gzip files
	found, err := regexp.MatchString(`\.[gG][zZ]$`, file)
	if err != nil {
		return fasta, fmt.Errorf("error matching gzip file pattern %w", err)
	}
	if found {
		// For gzip files, put a gzip.Reader into the chain
		gzr, err := gzip.NewReader(ff)
		if err != nil {
			return fasta, fmt.Errorf("unable to open gzip file %v: %w", file, err)
		}
		fasta.scanner = bufio.NewScanner(gzr)
	} else {
		// For non gzip files, go straight to bufio.Reader
		fasta.scanner = bufio.NewScanner(ff)
	}

	// Unnecessary but explicit
	fasta.scanner.Split(bufio.ScanLines)

	// Read the file
	for fasta.scanner.Scan() {
		line := strings.TrimSuffix(fasta.scanner.Text(), "\n")

		// TO DO - skip empty lines

		if faHeaderRex.MatchString(line) {
			fasta.Headers = append(fasta.Headers, line)
		} else if faIdRex.MatchString(line) {
			// Found Id line of the first record so save and return
			fasta.nextRecId = line
			return fasta, nil
		} else {
			return fasta, fmt.Errorf("should be impossible to get here - problematic line: %s", line)
		}
	}

	return fasta, nil
}

// Next returns the next record from the FASTA file. If there are no
// more records, it returns nil.
func (f *FastaFile) Next() (*FastaRec, error) {
	if f.EOF {
		return nil, nil
	}

	thisRec := NewFastaRec(f.nextRecId)
	thisRec.FastaFile = f
	f.recCtr++
	var seq strings.Builder

	for f.scanner.Scan() {
		line := strings.TrimSuffix(f.scanner.Text(), "\n")
		if faIdRex.MatchString(line) {
			f.nextRecId = line
			thisRec.Sequence = seq.String()
			return thisRec, nil
		} else {
			_, err := seq.WriteString(line)
			if err != nil {
				return thisRec, fmt.Errorf("error building sequence string: %w", err)
			}
		}
	}

	// Reached end-of-file
	thisRec.Sequence = seq.String()
	f.EOF = true
	return thisRec, nil
}

// ReadAll returns all of the remaining records from the FASTA file. If
// Next() has not been called then it will return all of the records
// but if Next() has been called, it will return the remaining records.
// If there are no records, it returns nil.
func (f *FastaFile) ReadAll() ([]*FastaRec, error) {
	var seqs []*FastaRec
	for {
		fr, err := f.Next()
		if err != nil {
			return seqs, err
		}
		if fr == nil {
			break
		}
        seqs = append(seqs,fr)
	}
	return seqs, nil
}

// RecordCount returns the number of records returned with Next().
func (f *FastaFile) RecordCount() int {
	return f.recCtr
}

// MD5 will return the MD5 string for the file and will calculate it on
// the first call, which is therefore slow. Subsequent calls return the
// already-calculated value.
func (f *FastaFile) MD5() (string, error) {
    //fmt.Printf("FastaFile:  %+v\n",f)
	if f.md5 != "" {
		return f.md5, nil
	}
	md5, err := Md5sum(f.Filepath)
	if err != nil {
		return "", fmt.Errorf("error generating MD5 for %v: %w", f.Filepath, err)
	}
	f.md5 = md5
	return f.md5, nil
}

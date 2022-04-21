package genome

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"regexp"
	"strings"

	log "github.com/sirupsen/logrus"
)

// FastaFile must be uniquely identifiable so we can check that a given
// FASTA file matches the one used to create Genomes etc. Files can
// change location so Filepath is not enough - MD5 is not a perfect
// solution but it's a pretty handy siganture.
type FastaFile struct {
	Filepath string
	MD5      string
}

// ParseFastaFile reads a file and parses it as though it were a FASTA
// file. It will handle gzipped files as long as they have a .gz
// extension.
func ParseFastaFile(file string) ([]*Sequence, error) {
	var seqs []*Sequence

	// Every FASTA file must be MD5 summed as part of the reading
	// process to provide a "signature" that we can use in future to
	// match FASTA files against serialised derived structures such as
	// Genome and Seed.
	md5, err := Md5sum(file)
	if err != nil {
		return nil, err
	}

	fasta := &FastaFile{
		Filepath: file,
		MD5:      md5}

	// Open file
	ff, err := os.Open(file)
	if err != nil {
		return seqs, err
	}
	defer ff.Close()

	// We need to define this before we handle gzip
	var scanner *bufio.Scanner

	// Based on file extension, handle gzip files
	found, err := regexp.MatchString(`\.[gG][zZ]$`, file)
	if err != nil {
		return seqs, fmt.Errorf("error matching gzip file pattern against %s: %w", file, err)
	}
	if found {
		// For gzip files, put a gzip.Reader into the chain
		reader, err := gzip.NewReader(ff)
		if err != nil {
			return seqs, fmt.Errorf("error opening gzip file %s: %w", file, err)
		}
		defer reader.Close()
		scanner = bufio.NewScanner(reader)
	} else {
		// For non gzip files, go straight to bufio.Reader
		scanner = bufio.NewScanner(ff)
	}

	// Unnecessary but explicit
	scanner.Split(bufio.ScanLines)

	// Pattern for header lines
	rex := regexp.MustCompile(`^>(.*)$`)

	// Read the file
	var thisSeq *Sequence
	var seqLines []string
	for scanner.Scan() {
		line := strings.TrimSuffix(scanner.Text(), "\n")
		if rex.MatchString(line) {
			// If not the first sequence, save to seqs
			if len(seqs) != 0 {
				// Use strings pkg to efficiently concatenate (potentially)
				// millions of lines (strings) containing sequence bases
				var builder strings.Builder
				for _, s := range seqLines {
					_, err := builder.WriteString(s)
					if err != nil {
						return seqs, fmt.Errorf("error building sequence string: %w", err)
					}
				}
				thisSeq.Sequence = builder.String()
				log.Debugf("  found sequence: (%d) %s", thisSeq.Length(), thisSeq.Header)
			}
			// (re)Initialise sequence reading machinery
			thisSeq = NewSequence(line)
			thisSeq.FastaFile = fasta
			seqs = append(seqs, thisSeq)
			seqLines = []string{}
		} else {
			// Defer building sequence string because, with immutable strings,
			// every concatenation is a memory allocation and copy - expensive!
			seqLines = append(seqLines, line)
		}
	}

	// Make sure we capture the last sequence
	var builder strings.Builder
	for _, s := range seqLines {
		_, err := builder.WriteString(s)
		if err != nil {
			return seqs, fmt.Errorf("error building sequence string: %w", err)
		}
	}
	thisSeq.Sequence = builder.String()
	log.Debugf("  found sequence: (%d) %s", thisSeq.Length(), thisSeq.Header)

	return seqs, nil
}

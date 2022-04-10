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

func ParseFastaFile(file string) ([]*Sequence, error) {
	var seqs []*Sequence

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
			line = strings.TrimPrefix(line, ">")
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
				log.Infof("  found sequence: %s %d", thisSeq.Header, thisSeq.Length())
			}
			// (re)Initialise sequence reading machinery
			thisSeq = NewSequence(line)
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

	return seqs, nil
}

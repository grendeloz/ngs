package genome

import (
	"bufio"
	"crypto/md5"
	"fmt"
	"io"
	"os"
	"strings"
)

// md5sum returns the MD5 hash of a file.  The MD5 provides a signature
// the file, allowing us to check whether two versions are the same.
func md5sum(file string) (string, error) {
	var chk string
	f, err := os.Open(file)
	if err != nil {
		return chk, err
	}
	defer f.Close()

	h := md5.New()
	if _, err := io.Copy(h, f); err != nil {
		return chk, err
	}

	chk = fmt.Sprintf("%x", h.Sum(nil))
	return chk, nil
}

// LinesFromFile reads a file and returns the trimmed lines.
func LinesFromFile(file string) ([]string, error) {
	var lines []string

	// Open file
	ff, err := os.Open(file)
	if err != nil {
		return lines, err
	}
	defer ff.Close()

	scanner := bufio.NewScanner(ff)
	scanner.Split(bufio.ScanLines)

	for scanner.Scan() {
		line := strings.TrimSuffix(scanner.Text(), "\n")
		lines = append(lines, line)
	}

	return lines, nil
}

func reverseString(s string) string {
	rns := []rune(s) // convert to rune
	for i, j := 0, len(rns)-1; i < j; i, j = i+1, j-1 {
		rns[i], rns[j] = rns[j], rns[i]
	}
	return string(rns)
}

func reverseBytes(b []byte) []byte {
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = b[j], b[i]
	}
	return b
}

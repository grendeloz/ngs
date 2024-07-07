package vcf

import (
	"strconv"
	"strings"
)

const missing string = "."

// Record holds a VCF record.
type Record struct {
	OrigStr string // string as read from file
	Chrom   string
	Pos     int
	Id      string
	Ref     string
	Alt     string
	Qual    float64
	Filter  string
	Info    string
	Format  string
	Samples []string
}

func (r Record) String() string {
	var ss []string

	// Append missing value for empty fields
	ss = append(ss, checkMissing(r.Chrom))
	ss = append(ss, checkMissing(strconv.Itoa(r.Pos)))
	ss = append(ss, checkMissing(r.Id))
	ss = append(ss, checkMissing(r.Ref))
	ss = append(ss, checkMissing(r.Alt))
	ss = append(ss, checkMissing(strconv.FormatFloat(r.Qual, 'E', -1, 64)))
	ss = append(ss, checkMissing(r.Filter))
	ss = append(ss, checkMissing(r.Info))

	if len(r.Samples) > 0 {
		ss = append(ss, checkMissing(r.Format))
		for _, x := range r.Samples {
			ss = append(ss, checkMissing(x))
		}
	}

	return strings.Join(ss, "\t")
}

func checkMissing(v string) string {
	if v == "" {
		return missing
	} else {
		return v
	}
}

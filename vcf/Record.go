package vcf

import (
	"fmt"
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

	// We need to handle Qual and Pos differently because they will be
	// converted from text to int/float and so if they were missing, we
	// need to put missing (not 0) back when we write them out
	QualMissing bool
	PosMissing  bool
}

func (r Record) String() string {
	var ss []string

	// Be careful to append missing value for empty fields
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

func RecordFromString(line string) (*Record, error) {
	r := &Record{}

	fields := strings.Split(line, "\t")
	if len(fields) < 8 {
		return r, fmt.Errorf("record has fewer than 8 fields: %s", line)
	}

	// Parse mandatory 8 fields
	r.Chrom = fields[0]
	r.Id = fields[2]
	r.Ref = fields[3]
	r.Alt = fields[4]
	r.Filter = fields[6]
	r.Info = fields[7]

	// Convert POS and QUAL. Need to handle missing character!
	if fields[1] == missing {
		r.Pos = 0
		r.PosMissing = true
	} else {
		i, err := strconv.Atoi(fields[1])
		if err != nil {
			return r, fmt.Errorf("cannot parse POS [%s] to int", fields[1])
		}
		r.Pos = i
	}
	if fields[5] == missing {
		r.Qual = 0
		r.QualMissing = true
	} else {
		f, err := strconv.ParseFloat(fields[5], 64)
		if err != nil {
			return r, fmt.Errorf("cannot parse QUAL [%s] to float", fields[5])
		}
		r.Qual = f
	}

	// Work out if this line contains genotypes in which case there
	// should also be a FORMAT field and one or more samples.

	if len(fields) > 8 {
		r.Format = fields[8]
		r.Samples = fields[9:]
	}

	return r, nil
}

func checkMissing(v string) string {
	if v == "" {
		return missing
	} else {
		return v
	}
}

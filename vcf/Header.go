package vcf

import (
	"strings"
)

type Header struct {
	OrigStr string   // string as read from file
	Samples []string // a VCF may have 0 samples
}

func NewHeader() *Header {
	return &Header{Samples: make([]string, 0, 4)}
}

func (h Header) String() string {
	// There are 8 mandatory columns plus FORMAT which only appears if
	// there are samples
	s := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
	if len(h.Samples) > 0 {
		s = s + "\tFORMAT\t"
		s = s + strings.Join(h.Samples, "\t")
	}

	return s
}

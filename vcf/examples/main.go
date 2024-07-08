package main

import (
	"fmt"
	"github.com/grendeloz/ngs/vcf"
)

func main() {
	stem := "/Users/johnp/DataTest/colo829/analysis_fbe3b136-dc8b-4c8d-bde3-a6390c91b521/"
	//file := stem + "jp.vcf.gz"
	file := stem + "fbe3b136-dc8b-4c8d-bde3-a6390c91b521.vcf.gz"
	v, err := vcf.NewFromFile(file)
	if err != nil {
		panic(err)
	}
	fmt.Printf("VCF file\n%+v", v)
	for i, r := range v.Records {
		if r.Chrom == `chr1` {
			fmt.Printf("%d:  %v  %v  %v %v %v\n", i, r.Chrom, r.Pos, r.Id, r.Ref, r.Alt)
		}
	}
}

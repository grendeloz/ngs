package genome

import (
	"testing"

	"github.com/google/go-cmp/cmp"
)

var read1 string = `@read1
ACTGGTAGTACTACTACATGGTCATTG
+read1
JJJBHDAGHJGGGHIJHFFFFDDHFFH
`

var read2 string = `@SRR357068.1 D042KACXX:3:1101:2690:2160 length=101
NCATCGTCCGGTATGTAGAACAGGGGAACCGGACGTTTTCCAAGGCGTAGCCATGTTAGACAAGGCGCAGATATAGGTGA
+SRR357068.1 D042KACXX:3:1101:2690:2160 length=101
#4=DBDDDHFHFFHIGHIIIJJJJJJJJJJJJBHDAGHJGGGHIJHFFFFDDEDCCDCCCCDDDDDBDBD>CDEE>C@CD
`

func TestNewRead(t *testing.T) {
	r1, err := NewReadFromString(read1)
	if err != nil {
		t.Fatalf(`NewReadFromString() failed: %v`, err)
	}

	r2, err := NewReadFromString(read2)
	if err != nil {
		t.Fatalf(`NewReadFromString() failed: %v`, err)
	}

	tests := []struct {
		name string
		want string
		got  string
	}{
		{name: "read1.Id", want: "read1", got: r1.Id},
		{name: "read1.Bases", want: "ACTGGTAGTACTACTACATGGTCATTG", got: string(r1.Bases)},
		{name: "read1.Qualities", want: "JJJBHDAGHJGGGHIJHFFFFDDHFFH", got: string(r1.Qualities)},
		{name: "read2.Id", want: "SRR357068.1 D042KACXX:3:1101:2690:2160 length=101", got: r2.Id},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			diff := cmp.Diff(tc.want, tc.got)
			if diff != "" {
				t.Fatalf(diff)
			}
		})
	}
}

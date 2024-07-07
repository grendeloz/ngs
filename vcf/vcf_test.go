package vcf

import (
	"strings"
	"testing"

	"github.com/google/go-cmp/cmp"
)

func TestVcf1(t *testing.T) {
	//file := ``
	tests := map[string]struct {
		input string
		sep   string
		want  []string
	}{
		"simple":       {input: "a/b/c", sep: "/", want: []string{"a", "b", "c"}},
		"wrong sep":    {input: "a/b/c", sep: ",", want: []string{"a/b/c"}},
		"no sep":       {input: "abc", sep: "/", want: []string{"abc"}},
		"trailing sep": {input: "a/b/c/", sep: "/", want: []string{"a", "b", "c", ""}},
	}

	for name, tc := range tests {
		t.Run(name, func(t *testing.T) {
			got := strings.Split(tc.input, tc.sep)
			diff := cmp.Diff(tc.want, got)
			if diff != "" {
				t.Fatalf(diff)
			}
		})
	}
}

package vcf

import (
	"github.com/grendeloz/kv"
)

// Holds Meta lines (anything that starts with '##').
type MetaRecord struct {
	IsStructured bool
	Key          string  // unstructured and structured lines
	Value        string  // unstructured lines only
	KVs          *kv.Set // structured lines only
}

func NewMetaRecord() *MetaRecord {
	return &MetaRecord{}
}

// Stringify Meta lines.
func (m *MetaRecord) String() string {
	if m.IsStructured {
		return `##` + m.Key + `=` + m.Value
	} else {
		return `##` + m.Key + `=<` + m.KVs.String() + `>`
	}
}

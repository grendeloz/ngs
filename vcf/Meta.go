package vcf

//import (
//    "fmt"
//)

// Holds Meta lines (anything that starts with '##').
type Meta struct {
	OrigStr string // string as read from file
	Records []*MetaRecord
}

func NewMeta() *Meta {
	return &Meta{Records: make([]*MetaRecord, 0, 10)}
}

// Stringify Meta lines.
func (m *Meta) String() string {
	return m.OrigStr
}

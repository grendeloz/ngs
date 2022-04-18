package genome

type Sequence struct {
	Header    string
	Sequence  string
	FastaFile *FastaFile
}

func NewSequence(header string) *Sequence {
	return &Sequence{Header: header}
}

func (s *Sequence) Length() int {
	return len(s.Sequence)
}

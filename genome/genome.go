package genome

import (
	"bufio"
	"encoding/gob"
	"os"

	"github.com/grendeloz/runp"
	"github.com/google/uuid"
	log "github.com/sirupsen/logrus"
)

// A Genome must be uniquely identifiable so we can check that derived
// objects, such as Seeds are only used with the correct Genome.
// This link between objects is especially important when Genomes are
// serialised to disk for later use.
type Genome struct {
	Name       string
	UUID       string
	Sequences  []*Sequence
	FastaFiles []*FastaFile
	Provenance []runp.RunParameters
}

func NewGenome(name string) *Genome {
	uuid := uuid.New()
	prov := []runp.RunParameters{runp.NewRunParameters()}
	return &Genome{
		Name:       name,
		UUID:       uuid.String(),
		Provenance: prov,
	}
}

type Sequence struct {
	Header    string
	Sequence  string
	FastaFile string
	FileMD5   string
}

func NewSequence(header string) *Sequence {
	return &Sequence{Header: header}
}

func (s *Sequence) Length() int {
    return len(s.Sequence)
}

// FastaFile must be uniquely identifiable so we can be sure in future
// that FASTA files match the FASTA files used to create Genomes,
// especially when Genomes have been serialised to disk.
type FastaFile struct {
	Name string
	MD5  string
}

// AddProvenance creates a new RunParameter and adds it onto the front
// (top) of the list of RunParameter in Provenance.
func (g *Genome) AddProvenance() {
	prov := runp.NewRunParameters()
	provs := []runp.RunParameters{prov}
	g.Provenance = append(provs, g.Provenance...)
}

func (g *Genome) NewSeed(seed string) (*Seed, error) {
	// Establish new Seed
	gs := &Seed{}
	gs.Mask = seed
	gs.genomeUUID = g.UUID
	gs.Offsets = make(map[string]int)
	gs.Coords = make(map[string][]int)

	// Set Provenance from source genome and then add new record
	gs.Provenance = g.Provenance
	gs.AddProvenance()

	// Add Provenance records from genome
	//for _, p := range g.Provenance {
	//	gs.Provenance = append(gs.Provenance, p)
	//}

	// Add Sequences from Genome
	for _, s := range g.Sequences {
		log.Infof("  adding sequence %s to Seed", s.Header)
		err := gs.addSequence(s)
		if err != nil {
			return gs, err
		}
	}

	// Apply Seed
	gs.applySeed(seed)

	return gs, nil
}

func (g *Genome) AddFastaFile(file string) error {
	// Every FASTA file must be MD5 summed as part of the reading
	// process to provide a "signature" that we can use in future to
	// match FASTA files against serialised derived structures such as
	// Genome and Seed.
	md5, err := Md5sum(file)
	if err != nil {
		return err
	}

	// Retrieve *Sequences from FASTA
	seqs, err := ParseFastaFile(file)
	if err != nil {
		return err
	}

	// Complete *Sequence information and add to Genome
	for _, s := range seqs {
		s.FastaFile = file
		s.FileMD5 = md5
		g.Sequences = append(g.Sequences, s)
	}

	// Add FASTA file to Genome
	ff := &FastaFile{
		Name: file,
		MD5:  md5}
	g.FastaFiles = append(g.FastaFiles, ff)

	return nil
}

// WriteAsGob serialises a genome to disk. The caller can specify the
// stem of the output filename but the pkg appends some identifying
// information.
func (g *Genome) WriteAsGob(filestem string) (string, error) {
	file := filestem + "." + g.UUID + ".genome.gob"

	f, err := os.Create(file)
	if err != nil {
		return file, err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	enc := gob.NewEncoder(w)
	err = enc.Encode(g)
	if err != nil {
		return file, err
	}
	return file, nil
}

func GenomeFromGob(file string) (*Genome, error) {
	var g *Genome

	// This is critical - gob will not decode to an empty (nil) pointer
	// type so we need to supply a real-but-empty variable.
	g = &Genome{}

	f, err := os.Open(file)
	if err != nil {
		return g, err
	}
	defer f.Close()

	r := bufio.NewReader(f)

	dec := gob.NewDecoder(r)
	err = dec.Decode(g)
	if err != nil {
		return g, err
	}
	return g, nil
}

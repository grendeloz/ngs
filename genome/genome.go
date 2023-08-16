// Package genome has types and functions for operating on genomics files
// and data.
package genome

import (
	"bufio"
	"encoding/gob"
	"fmt"
	"os"

	"github.com/google/uuid"
	"github.com/grendeloz/runp"
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
	Version    string
}

func NewGenome(name string) *Genome {
	uuid := uuid.New()
	prov := []runp.RunParameters{runp.NewRunParameters()}
	return &Genome{
		Name:       name,
		UUID:       uuid.String(),
		Provenance: prov,
		Version:    `0.1.0`,
	}
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
			return gs, fmt.Errorf("genome.Genome.NewSeed: %w", err)
		}
	}

	// Apply Seed
	gs.applySeed(seed)

	return gs, nil
}

func (g *Genome) AddFastaFile(file string) error {
	// Retrieve *Sequences from FASTA
	seqs, err := ParseFastaFile(file)
	if err != nil {
		return fmt.Errorf("genome.Genome.AddFastaFile: %w", err)
	}

	// Add to Genome
	g.Sequences = append(g.Sequences, seqs...)

	// Add FASTA file to Genome
	if len(seqs) > 0 {
		g.FastaFiles = append(g.FastaFiles, seqs[0].FastaFile)
	}

	return nil
}

// GetSequence returns a *Sequence or an error if the named sequence is
// not found. Note that the match is exact so case, spaces etc all
// matter - perfect match or no match.
func (g *Genome) GetSequence(seqName string) (*Sequence, error) {
	for _, s := range g.Sequences {
		if seqName == s.Name {
			return s, nil
		}
	}
	return nil, fmt.Errorf("Sequence %s not found in genome %s", seqName, g.Name)
}

// WriteAsGob serialises a genome to disk. The caller can specify the
// stem of the output filename but some identifying information is
// appended including the UUID. The filename is returned.
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

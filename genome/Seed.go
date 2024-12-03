package genome

import (
	"bufio"
	"encoding/gob"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/grendeloz/runp"
	log "github.com/sirupsen/logrus"
)

// Seed has a very different use case from Genome. It applies
// spaced seeds against a genome to create a list of locations where the
// seeded sequences appear.  A discussion of spaced seeds is beyond the
// scope of this documentation but the examples below should help.
// In essence a spaced seed is a way to create derived sequences by
// sampling non-contiguous positions from a source sequence. Only
// positions with a 1 in the seed are used to create the seed sequence.
// Spaced seeds are useful for finding inexact alignments without having
// a mismatch-tolerant algorithm. The seeds below show which bases from
// the source sequence would be kept in each seed and what the final
// seed sequence would be.
//
//	source  ACGTACGTACGTACGTACGT
//
//	1 seed  1_1_1_1_1_1_1_1_1_1_
//	        A G A G A G A G A G   =>  AGAGAGAGAG
//
//	2 seed  11__11__11__11__11__
//	        AC  AC  AC  AC  AC    =>  ACACACACAC
//
//	3 seed  ___1___1___1___1___1
//	           T   T   T   T   T  =>  TTTTT
//
//	4 seed  11111111111111111111
//	        ACGTACGTACGTACGTACGT  =>  ACGTACGTACGTACGTACGT
//
//	5 seed  1____1____1____1____
//	        A    C    G    T      =>  ACGT
//
// Because a Seed is fundamentally related to the Genome from
// which it is created, Seeds must be created via the Genome method
// NewSeed.
type Seed struct {
	Mask       string // e.g. 11_1_1
	Sequences  []*FastaRec
	Offsets    map[string]int
	Sequence   []byte
	Coords     map[string][]int
	Provenance []runp.RunParameters

	// This is intentionally private so it can only be accessed by
	// method GenomeUUID(). We don't want it to be user settable because we
	// want it to be an immutable record of the Genome that the
	// Seed came from.
	genomeUUID string
}

// For Sequence, []rune instead of []byte would obviously be
// preferable as it would let us cope with Unicode but the memory
// cost of rune (int32) vs byte (int8) is just not bearable for
// genome-scale (>3B records) slices.

// AddProvenance create a new RunParameter and adds it onto the front
// (top) of the list of RunParameter in Provenance.
func (gs *Seed) AddProvenance() {
	prov := runp.NewRunParameters()
	provs := []runp.RunParameters{prov}
	gs.Provenance = append(provs, gs.Provenance...)
}

// GenomeUUID returns the UUID assigned to this seed at creation.
func (gs *Seed) GenomeUUID() string {
	return gs.genomeUUID
}

// WriteAsGob serialises Seed in Go's gob binary format.
// The caller can set the output directory but cannot set the file name
// which has a fixed format. The name of the file written is returned.
func (gs *Seed) WriteAsGob(dir string) (string, error) {
	file := dir + "/" + gs.Mask + "." +
		gs.GenomeUUID() + ".seed.gob"

	f, err := os.Create(file)
	if err != nil {
		return file, err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	enc := gob.NewEncoder(w)
	err = enc.Encode(gs)
	if err != nil {
		return file, err
	}
	return file, nil
}

// SeedFromGob reads a file and unmarshals it assuming it to be a Seed
// serialised to disk using encoding/gob.
func SeedFromGob(file string) (*Seed, error) {
	var gs *Seed

	// This is critical - gob will not decode to an empty (nil) pointer
	// type so we need to supply a real-but-empty variable.
	gs = &Seed{}

	f, err := os.Open(file)
	if err != nil {
		return gs, err
	}
	defer f.Close()

	r := bufio.NewReader(f)

	dec := gob.NewDecoder(r)
	err = dec.Decode(gs)
	if err != nil {
		return gs, err
	}
	return gs, nil
}

// addSequence is a private function that only works to copy relevant
// pieces of a Sequence from a Genome to a Seed.
func (gs *Seed) addSequence(f *FastaRec) error {
	ns := Sequence{}

	// Note that Seed *Sequence do not hold the sequence bases -
	// these are in Seed.Sequence with Seed.Offsets
	// identifying where each Sequence starts and ends.
	ns.Header = f.Header
	ns.FastaFile = f.FastaFile

	// End of the current Seed sequence
	offset := len(gs.Sequence)

	gs.Offsets[s.Header] = offset
	gs.Sequences = append(gs.Sequences, &ns)
	gs.Sequence = append(gs.Sequence, []byte(f.Sequence)...)

	return nil
}

func (gs *Seed) applySeed(seed string) error {

	log.Infof("applying seed: %s", seed)
	// Work out which positions in the mask are interrogating.
	var seedpos []int
	seedlen := len(seed)
	seedbytes := []byte(seed)
	for i := 0; i < seedlen; i++ {
		if seedbytes[i] == '1' {
			seedpos = append(seedpos, i)
		}
	}
	log.Infof("  seed positions: %v", seedpos)
	seedposlen := len(seedpos)

	// To avoid lots of allocations we are going to use one seed slice
	// and overwrite it. Because strings are immutable, substring
	// replacements are pretty expensive so we will use a byte array
	// which we will rewrite and at the end we will string() it when we
	// need to use it as a hash lookup.
	var thisSeed []byte
	for i := 0; i < seedlen; i++ {
		thisSeed = append(thisSeed, '.')
	}

	// Apply the seed. For each sequence, construct the spaced seed at
	// every possible position and store the location in the uber-hash.
	lctr := 0
	for _, s := range gs.Sequences {
		log.Infof("  applying seed to: %s", s.Header)
		offset := gs.Offsets[s.Header]
		maxposn := offset + s.Length() - seedlen
		//log.Infof("    offset:%d  s.Length:%d  seedlen:%d maxposn:%d",
		//	offset, s.Length, seedlen, maxposn)
		for i := offset; i < maxposn; i++ {
			// This can be a progress reporter or a way to cut short
			// long chromosomes during testing
			if i%5000000 == 0 {
				log.Infof("    processing genomic position %d", i)
				break
			}
			//skip any oligo that starts with an N
			if gs.Sequence[i+seedpos[0]] == 'N' {
				continue
			}
			lctr++

			// Assemble oligo by walking the list of interrogating positions
			for j := 0; j < seedposlen; j++ {
				//fmt.Printf("%d %d %c %d %c\n", i, j, thisSeed[j], seedpos[j], mg.Sequence[i+seedpos[j]])
				thisSeed[j] = gs.Sequence[i+seedpos[j]]
			}
			// Add this position to the Offset - remember that we MUST
			// only take the first seedposlen bytes out of thisSeed which is
			// reused and so can be full of shite.

			oligo := string(thisSeed[0:seedposlen])
			gs.Coords[oligo] = append(gs.Coords[oligo], i)
		}
	}
	log.Infof("    Locations processed: %d", lctr)

	return nil
}

// WriteAsText serialises Seed as a text file. The text file is
// purely for debugging and does not attempt to write out all of the
// contents of a Seed.
// The caller can set the output directory but cannot set the file name
// which has a fixed format. The name of the file written is returned.
func (gs *Seed) WriteAsText(dir string) (string, error) {
	file := dir + "/" + gs.Mask + "." +
		gs.GenomeUUID() + ".seed.txt"

	f, err := os.Create(file)
	if err != nil {
		return file, err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	// Write Header
	maskheader := "# Seed: " + gs.Mask + "\n" +
		"# GenomeUUID: " + gs.genomeUUID + "\n"
	_, err = w.WriteString(maskheader)
	if err != nil {
		return file, fmt.Errorf("genome.Seed.WriteAsText: error writing header to %s: %w", maskheader, err)
	}

	// Write offsets - these are needed to interpret seed locations.
	// Let's go the extra mile and write them out in offset order.
	offsetMap := make(map[int]string)
	var offsets []int
	for s, i := range gs.Offsets {
		offsetMap[i] = s
		offsets = append(offsets, i)
	}
	sort.Ints(offsets)

	for _, i := range offsets {
		s := offsetMap[i]
		offset := fmt.Sprintf("# Offset,%s,%d\n", s, i)
		_, err := w.WriteString(offset)
		if err != nil {
			return file, fmt.Errorf("genome.Seed.WriteAsText: error writing offset %d: %w", i, err)
		}
	}

	// Write seeds and locations where they were found
	for seq, coords := range gs.Coords {
		var b strings.Builder

		// We know there is at least one coord so it simplifies the
		// separator handling if we manually handle the first coord and
		// then add any extras with separator chars.
		s := strconv.Itoa(coords[0])
		b.WriteString(seq + ":" + s)

		// Deal with any additional locations
		if len(coords) > 1 {
			for i := 1; i < len(coords); i++ {
				s := strconv.Itoa(coords[i])
				b.WriteString("," + s)
			}
		}

		// Write it all out
		_, err := w.WriteString(b.String() + "\n")
		if err != nil {
			return file, fmt.Errorf("genome.Seed.WriteAsText: error writing seed %s: %w", b.String(), err)
		}
	}

	return file, nil
}

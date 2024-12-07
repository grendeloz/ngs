# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

Changes since v0.4.0

## v0.4.0

### Changes
- Complete reorganisation of the genome package which introduced a very
large number of breaking changes - too many to sensibly list. You should
consider this release of genome to be a whole new package.
 - Handling of both FASTA and FASTQ files completely refactored to be
 more consistent and more efficient. 
 - Revised (and renamed) the record types for FASTA and FASTQ.
 - Sequence type is now a simple standalone type for sequence
 manipulation. Functions are provided to create Sequence from FastaRec
 and FastqRec types.

## v0.3.0

### Changes
- Changed to a single module at the root of the repo - the old model of
having multiple packages, each as a separate module, was awkward and was
causing needless trouble.

## v0.2.0

### Changes
- more detail in error messages to help with debugging.

### Additions
- Added simple Read type for next-gen sequencing reads.
- package gff3
- selector package
- Sequence.WithinLimits and test for same.

### Fixes
- updated logrus version to indirectly update a dependency with security
    problems.
- Fixed logic error in interval.OverlapsB case.
- Fixed bug in logic for creating 3 records from 2 original records for
    case where AOverlapsB.
- mistake in go.mod.

## v0.1.0

First release with the following major features:

genome package with code for
- reading FASTA files
- reading and writing genome objects in gob format
- early and incomplete code for spaced seeds

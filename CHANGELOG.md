# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

Changes since v0.2.0

### Changes
- Changed to a single module at the root of the repo - the old model of
having multiple packages, each as a separate module, was awkward and was
causing needless trouble.
- Complete reorganisation of the ngs package which introduced a very
large number of breaking changes. Too many to sensibly list - you should
consider this release of ngs to be a whole new package.
 - Handling of both FASTA and FASTQ files completely refactored to be
 more consistent and more efficient. 
 - Completely revised (and renamed) the record types for FASTA and
 FASTQ.
 - Sequence type is now a simple standalone type for sequence
 manipulation. Functions are provided to create Sequence from FastaRec
 and FastqRec types.

### Additions

### Fixes

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

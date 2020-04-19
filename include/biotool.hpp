// biotool
// Bioinformatics Toolbox
//
// include/biotool.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// API include file, provides access to biotool's features
//

#ifndef BIOTOOL_HPP
#define BIOTOOL_HPP

////////////////////////////////////////////////////////////////////////
// BEGIN FASTA

#include "../library/fasta.hpp"

#define FASTA Fasta
#define MOLECULE Fasta::Molecule
#define NUMBER_OF_MOLECULES(fasta) fasta.numberOfMolecules()
#define DESCRIPTION(molecule) molecule.getDescription()
#define SEQUENCE(molecule) molecule.getSequence()
#define SUBSEQUENCE(molecule, from, to) molecule.getSubSequence(from, to)
#define SEQUENCE_LENGTH(molecule) molecule.getSequenceLength()

// END FASTA
////////////////////////////////////////////////////////////////////////
// BEGIN PAIR ALIGNMENT

#include "../library/sequencePair.hpp"

#define SEQUENCE_PAIR SequencePair
#define SEQUENCES(sequencePair) sequencePair.getSequences()
#define HAMMING_DISTANCE(sequencePair) sequencePair.getHammingDistance()

// END PAIR ALIGNMENT
////////////////////////////////////////////////////////////////////////

#endif // BIOTOOL_HPP

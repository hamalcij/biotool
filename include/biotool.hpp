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
#define EDIT_DISTANCE(sequencePair) sequencePair.getEditDistance()

// END PAIR ALIGNMENT
////////////////////////////////////////////////////////////////////////
// BEGIN PDB

#include "../library/pdb.hpp"

#define PDB Pdb
#define MODEL Pdb::Model
#define CHAIN Pdb::Model::Chain
#define CHAINS std::vector<Pdb::Model::Chain>
#define RESIDUE Pdb::Model::Chain::Residue
#define RESIDUES std::vector<Pdb::Model::Chain::Residue>
#define HET_RESIDUE Pdb::Model::Chain::HetResidue
#define HET_RESIDUES std::vector<Pdb::Model::Chain::HetResidue>
#define ATOM Pdb::Model::Chain::Residue::Atom
#define ATOMS std::vector<Pdb::Model::Chain::Residue::Atom>
#define HET_ATOM Pdb::Model::Chain::HetResidue::HetAtom
#define HET_ATOMS std::vector<Pdb::Model::Chain::HetResidue::HetAtom>
#define GET_MODELS(pdb) pdb.getModels()
#define GET_CHAINS(model, chains) model.getChains(chains)
#define GET_RESIDUES(chain, residues) chain.getResidues(residues)
#define GET_HET_RESIDUES(chain, hetResidues) chain.getHetResidues(hetResidues)
#define GET_ATOMS(residue, atoms) residue.getAtoms(atoms)
#define FIND_MODEL(pdb, id) pdb.findModel(id)
#define FIND_CHAIN(model, id) model.findChain(id)
#define FIND_RESIDUE(chain, id) chain.findResidue(id)
#define FIND_HET_RESIDUE(chain, id) chain.findHetResidue(id)
#define FIND_ATOM(residue, id) residue.findAtom(id)
#define ID(unit) unit.getID()
#define RESIDUE_NAME(residue) residue.getResidueName()
#define ATOM_NAME(atom) atom.getAtomName()
#define COORDINATES(atom) atom.getCoords()
#define OCCUPANCY(atom) atom.getOcuupancy()
#define TEMPERATURE_FACTOR(atom) atom.getTempFactor()
#define ELEMENT(atom) atom.getElement()
#define CHARGE(atom) atom.getCharge()
#define NUMBER_OF_MODELS(pdb) pdb.getNumberOfModels()
#define NUMBER_OF_CHAINS(model) model.getNumberOfChains()
#define NUMBER_OF_RESIDUES(unit) unit.getNumberOfResidues()
#define NUMBER_OF_HET_RESIDUES(unit) unit.getNumberOfHetResidues()
#define NUMBER_OF_ATOMS(unit) unit.getNumberOfAtoms()
#define NUMBER_OF_HET_ATOMS(unit) unit.getNumberOfHetAtoms()

// END PDB
////////////////////////////////////////////////////////////////////////

#endif // BIOTOOL_HPP

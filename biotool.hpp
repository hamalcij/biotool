/*
MIT License

Copyright (c) 2020 Hamalčík Jan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// biotool
// Bioinformatics Toolbox
//
// biotool.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// API include file, provides access to biotool's features
//

#ifndef BIOTOOL_HPP
#define BIOTOOL_HPP

#include "library/fasta.hpp"
#include "library/sequencePair.hpp"
#include "library/pdb.hpp"

using namespace biotool;

////////////////////////////////////////////////////////////////////////
// BEGIN FASTA

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

#define SEQUENCE_PAIR SequencePair
#define SEQUENCES(sequencePair) sequencePair.getSequences()
#define HAMMING_DISTANCE(sequencePair) sequencePair.getHammingDistance()
#define EDIT_DISTANCE(sequencePair) sequencePair.getEditDistance()

// END PAIR ALIGNMENT
////////////////////////////////////////////////////////////////////////
// BEGIN PDB

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
#define GET_ATOMS_CLOSE_TO_LIGAND(model, ligand, atoms, maxDistance) model.getAtomsCloseToLigand(ligand, atoms, maxDistance)
#define GET_RESIDUES_CLOSE_TO_LIGAND(model, ligand, residues, maxDistance) model.getResiduesCloseToLigand(ligand, residues, maxDistance)
#define MODEL_WIDTH(model) model.getWidth()

// END PDB
////////////////////////////////////////////////////////////////////////

#endif // BIOTOOL_HPP

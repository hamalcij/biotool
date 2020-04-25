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
#include "library/clustal.hpp"

using namespace biotool;

////////////////////////////////////////////////////////////////////////
// BEGIN OVERLOADED

#define BT_CAT(A, B) A ## B
#define BT_EXPAND(...) __VA_ARGS__
#define BT_VA_ARGS_SIZE(...) BT_EXPAND(BT_APPLY_ARG_N((BT_ZERO_ARGS_DETECT(__VA_ARGS__), BT_RSEQ_N)))
#define BT_ZERO_ARGS_DETECT(...) BT_EXPAND(BT_ZERO_ARGS_DETECT_PREFIX_ ## __VA_ARGS__ ## _ZERO_ARGS_DETECT_SUFFIX)
#define BT_ZERO_ARGS_DETECT_PREFIX__ZERO_ARGS_DETECT_SUFFIX ,,,,,,,,,,,0
#define BT_APPLY_ARG_N(ARGS) BT_EXPAND(BT_ARGS_N ARGS)
#define BT_ARGS_N(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define BT_RSEQ_N 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
#define BT_OVERLOAD_SELECT(NAME, NUM) BT_CAT( NAME ## _, NUM)
#define BT_MACRO_OVERLOAD(NAME, ...) BT_OVERLOAD_SELECT(NAME, BT_VA_ARGS_SIZE(__VA_ARGS__))(__VA_ARGS__)

#define SEQUENCE(...) BT_MACRO_OVERLOAD(SEQUENCE, __VA_ARGS__)
#define ID(...) BT_MACRO_OVERLOAD(ID, __VA_ARGS__)
#define SUM_OF_PAIRS(...) BT_MACRO_OVERLOAD(SUM_OF_PAIRS, __VA_ARGS__)

// END OVERLOADED
////////////////////////////////////////////////////////////////////////
// BEGIN FASTA

#define FASTA Fasta
#define MOLECULE Fasta::Molecule
#define NUMBER_OF_MOLECULES(fasta) fasta.numberOfMolecules()
#define DESCRIPTION(molecule) molecule.getDescription()
#define SEQUENCE_1(molecule) molecule.getSequence()
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
#define ID_1(unit) unit.getID()
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
#define NUMBER_OF_SURFACE_AND_BURIED_RESIDUES(model) model.getNumberOfSurfaceAndBuried()
#define NUMBER_OF_HET_RESIDUES(unit) unit.getNumberOfHetResidues()
#define NUMBER_OF_ATOMS(unit) unit.getNumberOfAtoms()
#define NUMBER_OF_HET_ATOMS(unit) unit.getNumberOfHetAtoms()
#define GET_ATOMS_CLOSE_TO_LIGAND(model, ligand, atoms, maxDistance) model.getAtomsCloseToLigand(ligand, atoms, maxDistance)
#define GET_RESIDUES_CLOSE_TO_LIGAND(model, ligand, residues, maxDistance) model.getResiduesCloseToLigand(ligand, residues, maxDistance)
#define GET_SURFACE_AND_BURIED_STATS(model) model.getSurfaceAndBuriedStats()
#define MODEL_WIDTH(model) model.getWidth()
#define MODEL_DIAMETER(model) model.getDiameter()

// END PDB
////////////////////////////////////////////////////////////////////////
// BEGIN CLUSTAL

#define CLUSTAL Clustal
#define ID_2(clustal, index) clustal.getID(index)
#define SEQUENCE_2(clustal, specifier) clustal.getSequence(specifier)
#define COLUMN(clustal, index) clustal.getColumn(index)
#define NUMBER_OF_SEQUENCES(clustal) clustal.getNumberOfSequences()
#define NUMBER_OF_COLUMNS(clustal) clustal.getNumberOfColumns()
#define SUM_OF_PAIRS_3(clustal, matrix, order) clustal.sumOfPairs(matrix, order)
#define SUM_OF_PAIRS_4(clustal, matrix, order, index) clustal.sumOfPairs(matrix, order, index)
#define BEST_SCORING_COLUMNS(clustal, matrix, order, n) clustal.getBestScoringColumns(matrix, order, n)

// END CLUSTAL
////////////////////////////////////////////////////////////////////////

#endif // BIOTOOL_HPP

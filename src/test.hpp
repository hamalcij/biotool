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
// src/test.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// Test showing how biotool API works
//

#include "../biotool.hpp"

#include<iostream>

class TestFasta {
public:
  TestFasta(const std::string& file) : myFasta_(FASTA(file)) {}

  void performTest() {
    for (std::size_t i = 0; i < NUMBER_OF_MOLECULES(myFasta_); ++i) {
      MOLECULE myMolecule(myFasta_, i);
      std::cout << "NEW MOLECULE:" << std::endl;
      std::cout << "Description: " << DESCRIPTION(myMolecule) << std::endl;
      std::cout << "Sequence: " << SEQUENCE(myMolecule) << std::endl;
      std::cout << "Sequence length: "<< SEQUENCE_LENGTH(myMolecule) << std::endl;;
      std::cout << "Subsequence 2...24: " << SUBSEQUENCE(myMolecule, 2, 24) << std::endl;
      std::cout << std::endl;
    }
  }

private:
  FASTA myFasta_;
};

class TestSequencePair {
public:
  TestSequencePair(const std::string& seq1, const std::string& seq2) : mySequencePair_(SEQUENCE_PAIR(seq1, seq2)) {}

  void performTest() {
    auto[sequence1, sequence2] = SEQUENCES(mySequencePair_);
    auto[editDistance, alignment] = EDIT_DISTANCE(mySequencePair_);

    std::cout << "Sequence 1: " << sequence1 << std::endl;
    std::cout << "Sequence 2: " << sequence2 << std::endl;
    std::cout << "Hamming distance: " << HAMMING_DISTANCE(mySequencePair_) << std::endl;
    std::cout << "Edit distance: " << editDistance << std::endl;

    for (auto&& [seq1, seq2] : alignment) {
      std::cout << std::endl;
      std::cout << seq1 << std::endl;
      std::cout << seq2 << std::endl;
    }
    std::cout << std::endl;
  }

private:
  SEQUENCE_PAIR mySequencePair_;
};

class TestPdb {
public:
  TestPdb(const std::string& path) : myPdb_(PDB(path)) {}

  void performTest() {
    std::cout << "PDB" << std::endl;
    std::cout << "Models: " << NUMBER_OF_MODELS(myPdb_) << std::endl;

    for (auto&& model : GET_MODELS(myPdb_)) {
      std::cout << "Model " << ID(model);
      std::cout << "; Chains: " << NUMBER_OF_CHAINS(model);
      std::cout << "; Residues: " << NUMBER_OF_RESIDUES(model);
      std::cout << "; Atoms: " << NUMBER_OF_ATOMS(model);
      std::cout << "; HetResidues: " << NUMBER_OF_HET_RESIDUES(model);
      std::cout << "; HetAtoms: " << NUMBER_OF_HET_ATOMS(model);
      std::cout << std::endl;

      CHAINS chains;
      GET_CHAINS(model, chains);
      for (auto&& chain : chains) {
        std::cout << "  Chain " << ID(chain);
        std::cout << "; Residues: " << NUMBER_OF_RESIDUES(chain);
        std::cout << "; Atoms: " << NUMBER_OF_ATOMS(chain);
        std::cout << "; HetResidues: " << NUMBER_OF_HET_RESIDUES(chain);
        std::cout << "; HetAtoms: " << NUMBER_OF_HET_ATOMS(chain);
        std::cout << std::endl;

        RESIDUES residues;
        GET_RESIDUES(chain, residues);
        for (auto&& residue : residues) {
          std::cout << "    Residue " << RESIDUE_NAME(residue) << " " << ID(residue);
          std::cout << "; Atoms: " << NUMBER_OF_ATOMS(residue);
          std::cout << std::endl;

          ATOMS atoms;
          GET_ATOMS(residue, atoms);
          for (auto&& atom : atoms) {
            auto[x, y, z] = COORDINATES(atom);
            std::cout << "      Atom " << ATOM_NAME(atom) << " "<< ID(atom);
            std::cout << "; (" << x << ", " << y << ", " << z << ")";
            std::cout << "; Element: " << ELEMENT(atom);
            std::cout << "; Charge: " << CHARGE(atom);
            std::cout << std::endl;
          }
        }

        HET_RESIDUES hetResidues;
        GET_HET_RESIDUES(chain, hetResidues);
        for (auto&& hetResidue : hetResidues) {
          std::cout << "    HETResidue " << RESIDUE_NAME(hetResidue) << " " << ID(hetResidue);
          std::cout << "; Atoms: " << NUMBER_OF_HET_ATOMS(hetResidue);
          std::cout << std::endl;

          HET_ATOMS hetAtoms;
          GET_ATOMS(hetResidue, hetAtoms);
          for (auto&& hetAtom : hetAtoms) {
            auto[x, y, z] = COORDINATES(hetAtom);
            std::cout << "      HETAtom " << ATOM_NAME(hetAtom) << " "<< ID(hetAtom);
            std::cout << "; (" << x << ", " << y << ", " << z << ")";
            std::cout << "; Element: " << ELEMENT(hetAtom);
            std::cout << "; Charge: " << CHARGE(hetAtom);
            std::cout << std::endl;
          }
        }
      }
    }
    std::cout << std::endl;

    MODEL myModel = FIND_MODEL(myPdb_, "1");
    std::cout << "Found model " << ID(myModel);
    std::cout << "; Chains: " << NUMBER_OF_CHAINS(myModel);
    std::cout << "; Residues: " << NUMBER_OF_RESIDUES(myModel);
    std::cout << "; Atoms: " << NUMBER_OF_ATOMS(myModel);
    std::cout << "; HetResidues: " << NUMBER_OF_HET_RESIDUES(myModel);
    std::cout << "; HetAtoms: " << NUMBER_OF_HET_ATOMS(myModel);
    std::cout << std::endl;

    CHAIN myChain = FIND_CHAIN(myModel, 'B');
    std::cout << "  Found chain " << ID(myChain);
    std::cout << "; Residues: " << NUMBER_OF_RESIDUES(myChain);
    std::cout << "; Atoms: " << NUMBER_OF_ATOMS(myChain);
    std::cout << "; HetResidues: " << NUMBER_OF_HET_RESIDUES(myChain);
    std::cout << "; HetAtoms: " << NUMBER_OF_HET_ATOMS(myChain);
    std::cout << std::endl;

    RESIDUE myResidue = FIND_RESIDUE(myChain, "37");
    std::cout << "    Found residue " << RESIDUE_NAME(myResidue) << " " << ID(myResidue);
    std::cout << "; Atoms: " << NUMBER_OF_ATOMS(myResidue);
    std::cout << std::endl;

    ATOM myAtom = FIND_ATOM(myResidue, "1338");
    auto[x, y, z] = COORDINATES(myAtom);
    std::cout << "      Found atom " << ATOM_NAME(myAtom) << " " << ID(myAtom);
    std::cout << "; (" << x << ", " << y << ", " << z << ")";
    std::cout << "; Element: " << ELEMENT(myAtom);
    std::cout << "; Charge: " << CHARGE(myAtom);
    std::cout << std::endl;

    HET_RESIDUE myHetResidue = FIND_HET_RESIDUE(myChain, "147");
    std::cout << "    Found HET residue " << RESIDUE_NAME(myHetResidue) << " " << ID(myHetResidue);
    std::cout << "; Atoms: " << NUMBER_OF_HET_ATOMS(myHetResidue);
    std::cout << std::endl;

    ATOMS nearAtoms;
    float maxDistance = 3.2;
    GET_ATOMS_CLOSE_TO_LIGAND(myModel, myHetResidue, nearAtoms, maxDistance);
    for (auto&& atom : nearAtoms) {
      auto[nearX, nearY, nearZ] = COORDINATES(atom);
      std::cout << "      Near Atom within " << maxDistance << "Å ";
      std::cout << ATOM_NAME(atom) << " " << ID(atom);
      std::cout << "; Coords (" << nearX << ", " << nearY << ", " << nearZ << ")";
      std::cout << std::endl;
    }

    RESIDUES nearResidues;
    maxDistance = 2.6;
    GET_RESIDUES_CLOSE_TO_LIGAND(myModel, myHetResidue, nearResidues, maxDistance);
    for (auto&& residue : nearResidues) {
      std::cout << "      Near Residue within " << maxDistance << "Å ";
      std::cout << RESIDUE_NAME(residue) << " " << ID(residue);
      std::cout << "; Atoms: " << NUMBER_OF_ATOMS(residue);
      std::cout << std::endl;

      ATOMS nearResAtoms;
      GET_ATOMS(residue, nearResAtoms);
      for (auto&& atom : nearResAtoms) {
        auto[aX, aY, aZ] = COORDINATES(atom);
        std::cout << "        Atom " << ATOM_NAME(atom) << " " << ID(atom);
        std::cout << "; (" << aX << ", " << aY << ", " << aZ << ")";
        std::cout << "; Element: " << ELEMENT(atom);
        std::cout << "; Charge: " << CHARGE(atom);
        std::cout << std::endl;
      }
    }

    HET_ATOM myHetAtom = FIND_ATOM(myHetResidue, "4438");
    auto[hetX, hetY, hetZ] = COORDINATES(myHetAtom);
    std::cout << "      Found atom " << ATOM_NAME(myHetAtom) << " " << ID(myHetAtom);
    std::cout << "; (" << hetX << ", " << hetY << ", " << hetZ << ")";
    std::cout << "; Element: " << ELEMENT(myHetAtom);
    std::cout << "; Charge: " << CHARGE(myHetAtom);
    std::cout << std::endl;
  }

private:
  PDB myPdb_;
};

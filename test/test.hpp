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
// test/test.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// Test showing how biotool API works
//

#include "blosum64.hpp"
#include "../biotool.hpp"

#include<iostream>
#include<string>
#include<exception>

using namespace biotool;

class TestFasta {
public:
  TestFasta(const std::string& file) : myFasta_(FASTA(file)) {}

  void performTest() {
    int from, to;

    std::cout << std::endl;
    std::cout << "########################################################" << std::endl;
    std::cout << std::endl;
    std::cout << "TEST: FASTA" << std::endl;
    std::cout << "Type two sequence positions separated by space (e.g. 2 24). Subsequence spanning these positions will be fetched for each molecule:" << std::endl;
    std::cin >> from >> to;
    std::cout << std::endl;
    for (std::size_t i = 0; i < NUMBER_OF_MOLECULES(myFasta_); ++i) {
      MOLECULE myMolecule(myFasta_, i);
      std::cout << "MOLECULE:" << std::endl;
      std::cout << "Description: " << DESCRIPTION(myMolecule) << std::endl;
      std::cout << "Sequence: " << SEQUENCE(myMolecule) << std::endl;
      std::cout << "Sequence length: "<< SEQUENCE_LENGTH(myMolecule) << std::endl;;
      std::cout << "Subsequence " << from << "..." << to << ": " << SUBSEQUENCE(myMolecule, from, to) << std::endl;
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
    std::cout << std::endl;
    std::cout << "########################################################" << std::endl;
    std::cout << std::endl;
    std::cout << "TEST: Sequence pair" << std::endl;
    std::cout << std::endl;
    auto[sequence1, sequence2] = SEQUENCES(mySequencePair_);
    auto[editDistance, alignment] = EDIT_DISTANCE(mySequencePair_);

    std::cout << "Sequence 1: " << sequence1 << std::endl;
    std::cout << "Sequence 2: " << sequence2 << std::endl;
    try {
      std::cout << "Hamming distance: " << HAMMING_DISTANCE(mySequencePair_) << std::endl;
    }
    catch (...) {
      std::cout << "Cannot be computed, sequences do not have the same length!" << std::endl;
    }
    std::cout << "Edit distance: " << editDistance << std::endl;

    int i = 1;
    for (auto&& [seq1, seq2] : alignment) {
      std::cout << "Optimal alignment " << i << ":" << std::endl;
      std::cout << "  " << seq1 << std::endl;
      std::cout << "  " << seq2 << std::endl;
      ++i;
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
    std::string modelID, residueID, atomID, hetID, hetAtomID;
    char chainID;
    float maxDistance;

    std::cout << std::endl;
    std::cout << "########################################################" << std::endl;
    std::cout << std::endl;
    std::cout << "TEST: PDB" << std::endl;
    std::cout << "Type ID of a model (if the .pdb file doesn't specify a model, type simply \"1\"):" << std::endl;
    std::cin >> modelID;
    std::cout << "Type ID of a chain within model " << modelID << ":" << std::endl;
    std::cin >> chainID;
    std::cout << "Type ID of a residue within model " << modelID << " and chain " << chainID << ":" << std::endl;
    std::cin >> residueID;
    std::cout << "Type ID of an atom within model " << modelID << " and chain " << chainID << " and residue " << residueID << ":" << std::endl;
    std::cin >> atomID;
    std::cout << "Type ID of a ligand within model " << modelID << " and chain " << chainID << " (if there is no such ligand, just press enter):" << std::endl;
    std::cin >> hetID;
    if (!hetID.empty()) {
      std::cout << "Type ID of an atom within model " << modelID << " and chain " << chainID << " and ligand " << hetID << ":" << std::endl;
      std::cin >> hetAtomID;
      std::cout << "Type a floating point number representing a maximum distance to be searched for atoms and residues from ligand " << hetID << ":" << std::endl;
      std::cin >> maxDistance;
    }
    std::cout << std::endl;

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

    try {
      MODEL myModel = FIND_MODEL(myPdb_, modelID);
      std::cout << "Found model " << ID(myModel);
      std::cout << "; Chains: " << NUMBER_OF_CHAINS(myModel);
      std::cout << "; Residues: " << NUMBER_OF_RESIDUES(myModel);
      std::cout << "; Atoms: " << NUMBER_OF_ATOMS(myModel);
      std::cout << "; HetResidues: " << NUMBER_OF_HET_RESIDUES(myModel);
      std::cout << "; HetAtoms: " << NUMBER_OF_HET_ATOMS(myModel);
      std::cout << std::endl;

      try {
        CHAIN myChain = FIND_CHAIN(myModel, chainID);
        std::cout << "  Found chain " << ID(myChain);
        std::cout << "; Residues: " << NUMBER_OF_RESIDUES(myChain);
        std::cout << "; Atoms: " << NUMBER_OF_ATOMS(myChain);
        std::cout << "; HetResidues: " << NUMBER_OF_HET_RESIDUES(myChain);
        std::cout << "; HetAtoms: " << NUMBER_OF_HET_ATOMS(myChain);
        std::cout << std::endl;

        try {
          RESIDUE myResidue = FIND_RESIDUE(myChain, residueID);
          std::cout << "    Found residue " << RESIDUE_NAME(myResidue) << " " << ID(myResidue);
          std::cout << "; Atoms: " << NUMBER_OF_ATOMS(myResidue);
          std::cout << std::endl;

          try {
            ATOM myAtom = FIND_ATOM(myResidue, atomID);
            auto[x, y, z] = COORDINATES(myAtom);
            std::cout << "      Found atom " << ATOM_NAME(myAtom) << " " << ID(myAtom);
            std::cout << "; (" << x << ", " << y << ", " << z << ")";
            std::cout << "; Element: " << ELEMENT(myAtom);
            std::cout << "; Charge: " << CHARGE(myAtom);

            std::cout << std::endl;

            if (!hetID.empty()) {

              try {
                HET_RESIDUE myHetResidue = FIND_HET_RESIDUE(myChain, hetID);
                std::cout << "    Found HET residue " << RESIDUE_NAME(myHetResidue) << " " << ID(myHetResidue);
                std::cout << "; Atoms: " << NUMBER_OF_HET_ATOMS(myHetResidue);
                std::cout << std::endl;

                try {
                  ATOMS nearAtoms;
                  GET_ATOMS_CLOSE_TO_LIGAND(myModel, myHetResidue, nearAtoms, maxDistance);
                  for (auto&& atom : nearAtoms) {
                    auto[nearX, nearY, nearZ] = COORDINATES(atom);
                    std::cout << "      Near Atom within " << maxDistance << "Å ";
                    std::cout << ATOM_NAME(atom) << " " << ID(atom);
                    std::cout << "; Coords (" << nearX << ", " << nearY << ", " << nearZ << ")";
                    std::cout << std::endl;
                  }

                  try {
                    RESIDUES nearResidues;
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

                    try {
                      HET_ATOM myHetAtom = FIND_ATOM(myHetResidue, hetAtomID);
                      auto[hetX, hetY, hetZ] = COORDINATES(myHetAtom);
                      std::cout << "      Found atom " << ATOM_NAME(myHetAtom) << " " << ID(myHetAtom);
                      std::cout << "; (" << hetX << ", " << hetY << ", " << hetZ << ")";
                      std::cout << "; Element: " << ELEMENT(myHetAtom);
                      std::cout << "; Charge: " << CHARGE(myHetAtom);
                      std::cout << std::endl;
                    }
                    catch (const std::exception& e) {
                      std::cout << e.what() << std::endl;
                    }
                  }
                  catch (const std::exception& e) {
                    std::cout << e.what() << std::endl;
                  }
                }
                catch (const std::exception& e) {
                  std::cout << e.what() << std::endl;
                }
              }
              catch (const std::exception& e) {
                std::cout << e.what() << std::endl;
              }
            }
          }
          catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
          }
        }
        catch (const std::exception& e) {
          std::cout << e.what() << std::endl;
        }
      }
      catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
      }

      std::cout << "Model width: " << MODEL_WIDTH(myModel) << std::endl;
      std::cout << "Model diameter: " << MODEL_DIAMETER(myModel) << std::endl;
      std::cout << std::endl;

      auto [surface, buried] = NUMBER_OF_SURFACE_AND_BURIED_RESIDUES(myModel);
      std::cout << "Residues in model: " << NUMBER_OF_RESIDUES(myModel) << std::endl;
      std::cout << "Surface residues: " << surface << "; Buried residues: " << buried << std::endl;
      std::cout << std::endl;

      auto [surfaceStats, buriedStats] = GET_SURFACE_AND_BURIED_STATS(myModel);
      std::cout << "Statistics for surface residues:" << std::endl;
      for (auto&& [name, count] : surfaceStats) {
        std::cout << name << " " << count << std::endl;
      }
      std::cout << "Statistics for buried residues:" << std::endl;
      for (auto&& [name, count] : buriedStats) {
        std::cout << name << " " << count << std::endl;
      }

      auto [polarSurfacePortion, polarBuriedPortion] = GET_PORTION_OF_POLAR_SURFACE_AND_BURIED(myModel);
      std::cout << "Portion of polar residues on the surface: " << polarSurfacePortion;
      std::cout << "; Portion of polar residues in the core: " << polarBuriedPortion;
      std::cout << std::endl;
      std::cout << std::endl;
    }
    catch (const std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

private:
  PDB myPdb_;
};

class TestClustal {
public:
  TestClustal(const std::string& file) : myClustal_(CLUSTAL(file)) {}

  void performTest() {
    std::size_t index, n;

    std::cout << std::endl;
    std::cout << "########################################################" << std::endl;
    std::cout << std::endl;
    std::cout << "TEST: Clustal" << std::endl;
    std::cout << "Type index of a column in your MSA:" << std::endl;
    std::cin >> index;
    std::cout << "Type number of best scoring columns to be printed:" << std::endl;
    std::cin >> n;
    std::cout << std::endl;

    for (std::size_t i = 0; i < NUMBER_OF_SEQUENCES(myClustal_); ++i) {
      std::cout << ID(myClustal_, i) << ": " << SEQUENCE(myClustal_, i) << std::endl;
    }
    std::cout << std::endl;

    std::string myID = "UniRef90_UPI000";
    std::cout << "Found sequence of " << myID << std::endl;
    std::cout << SEQUENCE(myClustal_, myID) << std::endl;
    std::cout << std::endl;

    for (std::size_t i = 0; i < NUMBER_OF_COLUMNS(myClustal_); ++i) {
      for (auto&& aminoAcid : COLUMN(myClustal_, i)) {
        std::cout << aminoAcid;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Sum of pairs of column at index " << index << ": " << (int) SUM_OF_PAIRS(myClustal_, blosum64::triangularMatrix, blosum64::order, index) << std::endl;
    std::cout << "Sum of pairs with triangular matrix: " << (int) SUM_OF_PAIRS(myClustal_, blosum64::triangularMatrix, blosum64::order) << std::endl;
    std::cout << std::endl;

    std::cout << n << " best scoring columns:" << std::endl;
    for (auto&& [score, index, column] : BEST_SCORING_COLUMNS(myClustal_, blosum64::triangularMatrix, blosum64::order, n)) {
      std::cout << "Score: " << (int) score << "; Index: " << index;
      std::cout << "; Column: ";
      for (auto&& aminoAcid : *column) {
        std::cout << aminoAcid;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

private:
  CLUSTAL myClustal_;
};

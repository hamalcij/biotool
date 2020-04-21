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
// library/pdb.cpp
// Copyright (c) 2020 Hamalčík Jan
//
// Implementation of class Pdb
//

#include "pdb.hpp"

//
void Pdb::Model::Chain::Residue::getAtoms(Pdb::Model::Chain::Residue::atoms& atomVector) const {
  atomVector.clear();
  for (std::size_t i = indexFrom_; i <= indexTo_; ++i) {
    atomVector.emplace_back(model_, i);
  }
}

//
const Pdb::Model::Chain::Residue::Atom Pdb::Model::Chain::Residue::findAtom(const std::string& id) {
  for (std::size_t i = indexFrom_; i <= indexTo_; ++i) {
    if (model_.atomSerials_[i] == id) {
      return Pdb::Model::Chain::Residue::Atom(model_, i);
    }
  }
  throw std::invalid_argument("Atom ID " + id + " was not found!");
}

//
void Pdb::Model::Chain::HetResidue::getAtoms(Pdb::Model::Chain::HetResidue::hetAtoms& atomVector) const {
  atomVector.clear();
  for (std::size_t i = indexFrom_; i <= indexTo_; ++i) {
    atomVector.emplace_back(model_, i);
  }
}

//
const Pdb::Model::Chain::HetResidue::HetAtom Pdb::Model::Chain::HetResidue::findAtom(const std::string& id) {
  for (std::size_t i = indexFrom_; i <= indexTo_; ++i) {
    if (model_.hetAtomSerials_[i] == id) {
      return Pdb::Model::Chain::HetResidue::HetAtom(model_, i);
    }
  }
  throw std::invalid_argument("HET Atom ID " + id + " was not found!");
}

//
const std::size_t Pdb::Model::Chain::getNumberOfResidues() const {
  std::size_t numberOfResidues = 1;
  std::string current = model_.resSeqs_[indexFrom_];
  for (std::size_t i = indexFrom_ + 1; i < indexTo_; ++i) {
    if (current != model_.resSeqs_[i]) {
      current = model_.resSeqs_[i];
      ++numberOfResidues;
    }
  }

  return numberOfResidues;
}

//
const std::size_t Pdb::Model::Chain::getNumberOfHetResidues() const {
  if (model_.hetResSeqs_.empty()) {
    return 0;
  }
  std::set<std::string> uniqueResidues;
  auto resIt = model_.hetResSeqs_.cbegin();
  auto chainIt = model_.hetChainIDs_.cbegin();
  for (; resIt != model_.hetResSeqs_.cend(); ++resIt) {
    if (*chainIt == getID()) {
      uniqueResidues.emplace(*resIt);
    }
    ++chainIt;
  }

  return uniqueResidues.size();
}

//
const std::size_t Pdb::Model::Chain::getNumberOfHetAtoms() const {
  if (model_.hetChainIDs_.empty()) {
    return 0;
  }
  std::size_t numberOfHetAtoms = 0;
  for (auto&& chainID : model_.hetChainIDs_) {
    if (chainID == getID()) {
      ++numberOfHetAtoms;
    }
  }

  return numberOfHetAtoms;
}

//
void Pdb::Model::Chain::getResidues(Pdb::Model::Chain::residues& residueVector) const {
  std::size_t from = indexFrom_;
  std::string current = model_.resSeqs_[from];
  residueVector.clear();
  for (std::size_t i = from + 1; i <= indexTo_; ++i) {
    if (current != model_.resSeqs_[i]) {
      residueVector.emplace_back(model_, from, i - 1);
      from = i;
      current = model_.resSeqs_[i];
    }
  }
  residueVector.emplace_back(model_, from, indexTo_);
}

//
void Pdb::Model::Chain::getHetResidues(Pdb::Model::Chain::hetResidues& residueVector) const {
  std::size_t from = 0;
  char currentChain = model_.hetChainIDs_[from];
  while (from + 1 < model_.hetChainIDs_.size() && currentChain != getID()) {
    ++from;
    currentChain = model_.hetChainIDs_[from];
  }
  if (from == model_.hetChainIDs_.size()) {
    std::invalid_argument("Chain " + std::string(1, getID()) + " does not contain any HET residue!");
  }
  else {
    std::string current = model_.hetResSeqs_[from];
    bool isCorrectChain = true;
    residueVector.clear();
    for (std::size_t i = from + 1; i < model_.hetResSeqs_.size(); ++i) {
      if (model_.hetChainIDs_[i] == getID()) {
        if (isCorrectChain) {
          if (current != model_.hetResSeqs_[i]) {
            residueVector.emplace_back(model_, from, i - 1);
            from = i;
            current = model_.hetResSeqs_[i];
          }
        }
        else {
          isCorrectChain = true;
          from = i;
          current = model_.hetResSeqs_[i];
        }
      }
      else {
        if (isCorrectChain) {
          residueVector.emplace_back(model_, from, i - 1);
          isCorrectChain = false;
        }
      }
    }
    if (isCorrectChain) {
      residueVector.emplace_back(model_, from, model_.hetResSeqs_.size() - 1);
    }
  }
}

//
const Pdb::Model::Chain::Residue Pdb::Model::Chain::findResidue(const std::string& id) const {
  std::size_t from = indexFrom_;
  bool notFound = true;
  for (std::size_t i = from; i <= indexTo_; ++i) {
    if (model_.resSeqs_[i] == id) {
      if (notFound) {
        from = i;
        notFound = false;
      }
    }
    else if (!notFound) {
      return Pdb::Model::Chain::Residue(model_, from, i - 1);
    }
  }
  if (notFound) {
    throw std::invalid_argument("Residue ID " + id + " was not found!");
  }
  else {
    return Pdb::Model::Chain::Residue(model_, from, indexTo_);
  }
}

//
const Pdb::Model::Chain::HetResidue Pdb::Model::Chain::findHetResidue(const std::string& id) const {
  if (model_.hetResSeqs_.empty()) {
    throw std::invalid_argument("Chain " + std::string(1, getID()) + " does not contain any HET residue!");
  }
  std::size_t from = 0;
  bool notFound = true;
  for (std::size_t i = 0; i < model_.hetResSeqs_.size(); ++i) {
    if (model_.hetChainIDs_[i] == getID() && model_.hetResSeqs_[i] == id) {
      if (notFound) {
        from = i;
        notFound = false;
      }
    }
    else if (!notFound) {
      return Pdb::Model::Chain::HetResidue(model_, from, i - 1);
    }
  }
  if (notFound) {
    throw std::invalid_argument("HET Residue ID " + id + " was not found!");
  }
  else {
    return Pdb::Model::Chain::HetResidue(model_, from, model_.hetResSeqs_.size() - 1);
  }
}

//
void Pdb::Model::getAtomsCloseToLigand(const Pdb::Model::Chain::HetResidue& het, Pdb::Model::Chain::Residue::atoms& atoms, const float maxDistance) const {
  Pdb::Model::Chain::HetResidue::hetAtoms hetAtoms;
  het.getAtoms(hetAtoms);
  atoms.clear();

  for (std::size_t i = 0; i < coords_.size(); ++i) {
    auto[x, y, z] = coords_[i];
    for (auto&& hetAtom : hetAtoms) {
      auto[hetX, hetY, hetZ] = hetAtom.getCoords();
      if (distance(x, y, z, hetX, hetY, hetZ) <= maxDistance) {
        atoms.emplace_back(*this, i);
        break;
      }
    }
  }
}

//
void Pdb::Model::getResiduesCloseToLigand(const Pdb::Model::Chain::HetResidue& het, Pdb::Model::Chain::residues& residues, const float maxDistance) const {
  Pdb::Model::Chain::HetResidue::hetAtoms hetAtoms;
  het.getAtoms(hetAtoms);
  residues.clear();

  Pdb::Model::chains chains;
  getChains(chains);
  for (auto&& chain : chains) {

    Pdb::Model::Chain::residues chainResidues;
    chain.getResidues(chainResidues);
    for (auto&& residue : chainResidues) {
      bool residueIsNear = false;

      Pdb::Model::Chain::Residue::atoms atoms;
      residue.getAtoms(atoms);
      for (auto&& atom : atoms) {

        auto[x, y, z] = atom.getCoords();
        for (auto&& hetAtom : hetAtoms) {
          auto[hetX, hetY, hetZ] = hetAtom.getCoords();
          if (distance(x, y, z, hetX, hetY, hetZ) <= maxDistance) {
            residues.push_back(residue);
            residueIsNear = true;
            break;
          }
        }
        if (residueIsNear) break;
      }
    }
  }
}

//
void Pdb::Model::getChains(Pdb::Model::chains& chainVector) const {
  std::size_t from = 0;
  char current = chainIDs_[from];
  chainVector.clear();
  for (std::size_t i = from + 1; i < chainIDs_.size(); ++i) {
    if (current != chainIDs_[i]) {
      chainVector.emplace_back(*this, from, i - 1);
      from = i;
      current = chainIDs_[i];
    }
  }
  chainVector.emplace_back(*this, from, chainIDs_.size() - 1);
}

//
const Pdb::Model::Chain Pdb::Model::findChain(const char id) const {
  std::size_t from = 0;
  bool notFound = true;
  for (std::size_t i = 0; i < chainIDs_.size(); ++i) {
    if (chainIDs_[i] == id) {
      if (notFound) {
        from = i;
        notFound = false;
      }
    }
    else if (!notFound) {
      return Pdb::Model::Chain(*this, from, i - 1);
    }
  }
  if (notFound) {
    throw std::invalid_argument("Chain ID " + std::string(1, id) + " was not found!");
  }
  else {
    return Pdb::Model::Chain(*this, from, chainIDs_.size() - 1);
  }
}

//
const Pdb::Model& Pdb::findModel(const std::string& id) const {
  auto it = std::find_if(
    models_.cbegin(),
    models_.cend(),
    [id](const Pdb::Model& model){
      return id == model.id_;
    }
  );
  if (it == models_.cend()) {
    throw std::invalid_argument("Model ID " + id + " was not found!");
  }

  return *it;
}

//
void Pdb::parsePdbFile(const std::string& path) {
  std::ifstream file(path);

  if (file.good() && !fileIsEmpty(file)) {
    char line[bufferSize_];
    bool containsModels = false;
    bool awaitsENDMDL = false;
    bool awaitsAtomTER = false;
    bool awaitsHetTER = false;
    std::size_t currentModelIndex = 0;

    for (;;) {
      file.getline(line, bufferSize_);
      if (file.fail()) break;

      std::string section = getKeyword(line, recordNameFrom_, recordNameTo_);

      if (section == "MODEL") {
        if (awaitsENDMDL || awaitsAtomTER || awaitsHetTER) {
          pdbFileCorrupted(path);
        }
        else {
          containsModels = true;
          awaitsENDMDL = true;
          std::string modelID = getKeyword(line, modelSerialFrom_, modelSerialTo_);
          if (modelID.empty()) {
            pdbFileCorrupted(path);
          }
          models_.emplace_back(*this, modelID);
        }
      }

      else if (section == "ENDMDL") {
        if (awaitsENDMDL) {
          awaitsENDMDL = false;
          ++currentModelIndex;
        }
        else {
          pdbFileCorrupted(path);
        }
      }

      else if (section == "ATOM") {
        if (!awaitsHetTER) {
          if (!containsModels) {
            containsModels = true;
            models_.emplace_back(*this, "1");
          }
          awaitsAtomTER = true;
          try {
            models_[currentModelIndex].atomSerials_.emplace_back(
              getKeyword(line, atomSerialFrom_, atomSerialTo_)
            );
            models_[currentModelIndex].atomNames_.emplace_back(
              getKeyword(line, atomNameFrom_, atomNameTo_)
            );
            models_[currentModelIndex].altLocs_.emplace_back(
              line[altLocPos_]
            );
            models_[currentModelIndex].resNames_.emplace_back(
              getKeyword(line, resNameFrom_, resNameTo_)
            );
            models_[currentModelIndex].chainIDs_.emplace_back(
              line[chainIDPos_]
            );
            models_[currentModelIndex].resSeqs_.emplace_back(
              getKeyword(line, resSeqFrom_, resSeqTo_)
            );
            models_[currentModelIndex].iCodes_.emplace_back(
              line[iCodePos_]
            );
            models_[currentModelIndex].coords_.emplace_back(
              std::make_tuple(
                std::stof(getKeyword(line, xFrom_, xTo_)),
                std::stof(getKeyword(line, yFrom_, yTo_)),
                std::stof(getKeyword(line, zFrom_, zTo_))
              )
            );
            models_[currentModelIndex].occupancies_.emplace_back(
              std::stof(getKeyword(line, occupancyFrom_, occupancyTo_))
            );
            models_[currentModelIndex].tempFactors_.emplace_back(
              std::stof(getKeyword(line, tempFactorFrom_, tempFactorTo_))
            );
            models_[currentModelIndex].elements_.emplace_back(
              getKeyword(line, elementFrom_, elementTo_)
            );
            models_[currentModelIndex].charges_.emplace_back(
              getKeyword(line, chargeFrom_, chargeTo_)
            );
          }
          catch (...) {
            pdbFileCorrupted(path);
          }
        }
        else {
          pdbFileCorrupted(path);
        }
      }

      else if (section == "HETATM") {
        if (!awaitsAtomTER) {
          if (!containsModels) {
            containsModels = true;
            models_.emplace_back(*this, "");
          }
          awaitsHetTER = true;
          try {
            models_[currentModelIndex].hetAtomSerials_.emplace_back(
              getKeyword(line, atomSerialFrom_, atomSerialTo_)
            );
            models_[currentModelIndex].hetAtomNames_.emplace_back(
              getKeyword(line, atomNameFrom_, atomNameTo_)
            );
            models_[currentModelIndex].hetAltLocs_.emplace_back(
              line[altLocPos_]
            );
            models_[currentModelIndex].hetResNames_.emplace_back(
              getKeyword(line, resNameFrom_, resNameTo_)
            );
            models_[currentModelIndex].hetChainIDs_.emplace_back(
              line[chainIDPos_]
            );
            models_[currentModelIndex].hetResSeqs_.emplace_back(
              getKeyword(line, resSeqFrom_, resSeqTo_)
            );
            models_[currentModelIndex].hetICodes_.emplace_back(
              line[iCodePos_]
            );
            models_[currentModelIndex].hetCoords_.emplace_back(
              std::make_tuple(
                std::stof(getKeyword(line, xFrom_, xTo_)),
                std::stof(getKeyword(line, yFrom_, yTo_)),
                std::stof(getKeyword(line, zFrom_, zTo_))
              )
            );
            models_[currentModelIndex].hetOccupancies_.emplace_back(
              std::stof(getKeyword(line, occupancyFrom_, occupancyTo_))
            );
            models_[currentModelIndex].hetTempFactors_.emplace_back(
              std::stof(getKeyword(line, tempFactorFrom_, tempFactorTo_))
            );
            models_[currentModelIndex].hetElements_.emplace_back(
              getKeyword(line, elementFrom_, elementTo_)
            );
            models_[currentModelIndex].hetCharges_.emplace_back(
              getKeyword(line, chargeFrom_, chargeTo_)
            );
          }
          catch (...) {
            pdbFileCorrupted(path);
          }
        }
        else {
          pdbFileCorrupted(path);
        }
      }

      else if (section == "TER") {
        if (awaitsAtomTER) {
          awaitsAtomTER = false;
        }
        else if (awaitsHetTER) {
          awaitsHetTER = false;
        }
      }
    }
  }
  else {
    pdbFileCorrupted(path);
  }

  file.close();
}

//
std::string Pdb::getKeyword(
  char* line,
  const unsigned short from,
  const unsigned short to
) {
  std::string keyword = "";

  for (unsigned short i = from; i <= to; ++i) {
    if (!isspace(line[i])) {
      keyword += line[i];
    }
  }

  return keyword;
}

//
void Pdb::pdbFileCorrupted(const std::string& path) {
  throw std::ifstream::failure("PDB file at " + path + " is corrupted!");
}

//
bool Pdb::fileIsEmpty(std::ifstream& file) {
  return file.peek() == std::ifstream::traits_type::eof();
}

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

namespace biotool {

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
  const Pdb::Model::size_tTuple Pdb::Model::getNumberOfSurfaceAndBuried() {
    getSurfaceResidues();

    return std::make_tuple(
      surfaceResidues_.size(),
      getNumberOfResidues() - surfaceResidues_.size()
    );
  }

  //
  const Pdb::Model::statsTuple Pdb::Model::getSurfaceAndBuriedStats() {
    if (surfaceStats_.empty()) {
      getSurfaceResidues();

      Pdb::Model::residuesSet buried;
      Pdb::Model::chains chains;
      getChains(chains);
      for (auto&& chain : chains) {
        Pdb::Model::Chain::residues residues;
        chain.getResidues(residues);
        for (auto&& residue : residues) {
          auto search = surfaceResidues_.find(residue);
          if (search == surfaceResidues_.end()) {
            buried.insert(residue);
          }
        }
      }

      surfaceStats_.insert({"ALA", 0});
      surfaceStats_.insert({"ARG", 0});
      surfaceStats_.insert({"ASN", 0});
      surfaceStats_.insert({"ASP", 0});
      surfaceStats_.insert({"CYS", 0});
      surfaceStats_.insert({"GLN", 0});
      surfaceStats_.insert({"GLU", 0});
      surfaceStats_.insert({"GLY", 0});
      surfaceStats_.insert({"HIS", 0});
      surfaceStats_.insert({"ILE", 0});
      surfaceStats_.insert({"LEU", 0});
      surfaceStats_.insert({"LYS", 0});
      surfaceStats_.insert({"MET", 0});
      surfaceStats_.insert({"PHE", 0});
      surfaceStats_.insert({"PRO", 0});
      surfaceStats_.insert({"SER", 0});
      surfaceStats_.insert({"THR", 0});
      surfaceStats_.insert({"TRP", 0});
      surfaceStats_.insert({"TYR", 0});
      surfaceStats_.insert({"VAL", 0});
      buriedStats_ = surfaceStats_;

      for (auto&& residue : surfaceResidues_) {
        ++surfaceStats_[residue.getResidueName()];
      }
      for (auto&& residue : buried) {
        ++buriedStats_[residue.getResidueName()];
      }
    }

    return std::make_tuple(surfaceStats_, buriedStats_);
  }

  //
  const Pdb::Model::floatTuple Pdb::Model::getPortionOfPolarResidues() {
    std::set<std::string> polarResidues = {
      "ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "LYS", "SER", "THR", "TYR"
    };

    getSurfaceAndBuriedStats();
    auto [surfaceCount, buriedCount] = getNumberOfSurfaceAndBuried();
    std::size_t polarSurfaceCount = 0;
    std::size_t polarBuriedCount = 0;

    for (auto&& polar : polarResidues) {
      polarSurfaceCount += surfaceStats_[polar];
      polarBuriedCount += buriedStats_[polar];
    }

    return std::make_tuple(
      static_cast<float>(polarSurfaceCount) / static_cast<float>(surfaceCount),
      static_cast<float>(polarBuriedCount) / static_cast<float>(buriedCount)
    );
  }

  //
  const float Pdb::Model::getWidth() {
    if (std::get<0>(farthestAtoms_) == 0) {
      getFarthestAtoms();
    }

    return std::get<0>(farthestAtoms_);
  }

  //
  const float Pdb::Model::getDiameter() {
    if (diameter_ == 0) {
      getFarthestAtoms();
      auto&[diam, atom1, atom2] = farthestAtoms_;
      float radius = diam / 2.0;
      Pdb::Model::fVector3 center(
        (atom1.x + atom2.x) / 2.0,
        (atom1.y + atom2.y) / 2.0,
        (atom1.z + atom2.z) / 2.0
      );

      for (;;) {
        Pdb::Model::fVector3 mostDistantAtom;
        float maxDistance = 0;
        for (auto&& atom : convexHull_) {
          const float dist = atom.getDistanceTo(center);
          if (dist > maxDistance) {
            maxDistance = dist;
            mostDistantAtom.x = atom.x;
            mostDistantAtom.y = atom.y;
            mostDistantAtom.z = atom.z;
          }
        }
        if (maxDistance > radius) {
          center = circumcenter(atom1, atom2, mostDistantAtom);
          radius *= 1.00001;
        }
        else {
          diameter_ = radius * 2.0;
          break;
        }
      }
    }

    return diameter_;
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
  Pdb::Model::fVector3 Pdb::Model::circumcenter(
    const Pdb::Model::fVector3& a,
    const Pdb::Model::fVector3& b,
    const Pdb::Model::fVector3& c
  ) {

    const float ab = a.getDistanceTo(b);
    const float bc = b.getDistanceTo(c);
    const float ca = c.getDistanceTo(a);

    const float alpha = bc*bc * (ca*ca + ab*ab - bc*bc);
    const float beta = ca*ca * (ab*ab + bc*bc - ca*ca);
    const float gamma = ab*ab * (bc*bc + ca*ca - ab*ab);

    const float x = (alpha*a.x + beta*b.x + gamma*c.x) / (alpha + beta + gamma);
    const float y = (alpha*a.y + beta*b.y + gamma*c.y) / (alpha + beta + gamma);
    const float z = (alpha*a.z + beta*b.z + gamma*c.z) / (alpha + beta + gamma);

    return Pdb::Model::fVector3(x, y, z);
  }

  //
  void Pdb::Model::getFarthestAtoms() {
    if (createConvexHull()) {
      for (std::size_t i = 0; i < convexHull_.size(); ++i) {
        for (std::size_t j = i + 1; j < convexHull_.size(); ++j) {
          const float dist = convexHull_[i].getDistanceTo(convexHull_[j]);
          if (dist > std::get<0>(farthestAtoms_)) {
            std::get<0>(farthestAtoms_) = dist;
            std::get<1>(farthestAtoms_) = convexHull_[i];
            std::get<2>(farthestAtoms_) = convexHull_[j];
          }
        }
      }
    }
  }

  //
  const bool Pdb::Model::createConvexHull() {
    if (convexHull_.size() == 0) {
      using namespace quickhull;

      QuickHull<float> qh;
      std::vector<Vector3<float>> pointCloud;
      for (auto&& [x, y, z] : coords_) {
        pointCloud.emplace_back(x, y, z);
      }

      auto hull = qh.getConvexHull(pointCloud, false, false);
      convexHull_ = hull.getVertexBuffer();

      return true;
    }

    return false;
  }

  //
  const bool Pdb::Model::atomIsSideChain(const Chain::Residue::Atom& atom) {
    const std::string& atomName = atom.getAtomName();

    if (atomName == "N") return false;
    if (atomName == "CA") return false;
    if (atomName == "C") return false;
    if (atomName == "O") return false;
    if (atomName == "OXT") return false;
    return true;
  }

  //
  void Pdb::Model::getSurfaceResidues() {
    if (surfaceResidues_.empty()) {
      using cell = std::optional<std::tuple<Pdb::Model::Chain::Residue, Pdb::Model::Chain::Residue::Atom>>;
      using z = std::vector<cell>;
      using yz = std::vector<z>;
      using xyz = std::vector<yz>;
      using cellCoords = std::tuple<std::size_t, std::size_t, std::size_t>;

      const auto [minX, minY, minZ] = getMinCoords();
      const float cellSize = solvent_;
      const std::size_t gridSize = static_cast<std::size_t>(getDiameter() / cellSize) + 10;
      bool discovered[gridSize][gridSize][gridSize];

      xyz protein;
      for (std::size_t i = 0; i < gridSize; ++i) {
        yz yCoord;
        for (std::size_t j = 0; j < gridSize; ++j) {
          yCoord.emplace_back(gridSize);
        }
        protein.push_back(std::move(yCoord));
      }

      Pdb::Model::chains chains;
      getChains(chains);
      for (auto&& chain : chains) {
        Pdb::Model::Chain::residues residues;
        chain.getResidues(residues);
        for (auto&& residue : residues) {
          Pdb::Model::Chain::Residue::atoms atoms;
          residue.getAtoms(atoms);
          for (auto&& atom : atoms) {
            auto [x, y, z] = atom.getCoords();
            std::size_t xPos = static_cast<int>((x - minX) / cellSize) + 5;
            std::size_t yPos = static_cast<int>((y - minY) / cellSize) + 5;
            std::size_t zPos = static_cast<int>((z - minZ) / cellSize) + 5;
            protein[xPos][yPos][zPos] = std::make_tuple(residue, atom);
          }
        }
      }

      std::deque<cellCoords> queue;
      queue.emplace_front(std::make_tuple(0, 0, 0));
      discovered[0][0][0] = true;
      while (!queue.empty()) {
        auto [x, y, z] = queue.back();
        queue.pop_back();
        if (x > 0) {
          if (protein[x - 1][y][z].has_value()) {
            auto [residue, atom] = *protein[x - 1][y][z];
            if (atomIsSideChain(atom)) surfaceResidues_.insert(residue);
          }
          else if (!discovered[x - 1][y][z]) {
            discovered[x - 1][y][z] = true;
            queue.emplace_front(std::make_tuple(x - 1, y, z));
          }
        }
        if (x < gridSize - 1) {
          if (protein[x + 1][y][z].has_value()) {
            auto [residue, atom] = *protein[x + 1][y][z];
            if (atomIsSideChain(atom)) surfaceResidues_.insert(residue);
          }
          else if (!discovered[x + 1][y][z]) {
            discovered[x + 1][y][z] = true;
            queue.emplace_front(std::make_tuple(x + 1, y, z));
          }
        }
        if (y > 0) {
          if (protein[x][y - 1][z].has_value()) {
            auto [residue, atom] = *protein[x][y - 1][z];
            if (atomIsSideChain(atom)) surfaceResidues_.insert(residue);
          }
          else if (!discovered[x][y - 1][z]) {
            discovered[x][y - 1][z] = true;
            queue.emplace_front(std::make_tuple(x, y - 1, z));
          }
        }
        if (y < gridSize - 1) {
          if (protein[x][y + 1][z].has_value()) {
            auto [residue, atom] = *protein[x][y + 1][z];
            if (atomIsSideChain(atom)) surfaceResidues_.insert(residue);
          }
          else if (!discovered[x][y + 1][z]) {
            discovered[x][y + 1][z] = true;
            queue.emplace_front(std::make_tuple(x, y + 1, z));
          }
        }
        if (z > 0) {
          if (protein[x][y][z - 1].has_value()) {
            auto [residue, atom] = *protein[x][y][z - 1];
            if (atomIsSideChain(atom)) surfaceResidues_.insert(residue);
          }
          else if (!discovered[x][y][z - 1]) {
            discovered[x][y][z - 1] = true;
            queue.emplace_front(std::make_tuple(x, y, z - 1));
          }
        }
        if (z < gridSize - 1) {
          if (protein[x][y][z + 1].has_value()) {
            auto [residue, atom] = *protein[x][y][z + 1];
            if (atomIsSideChain(atom)) surfaceResidues_.insert(residue);
          }
          else if (!discovered[x][y][z + 1]) {
            discovered[x][y][z + 1] = true;
            queue.emplace_front(std::make_tuple(x, y, z + 1));
          }
        }
      }
    }
  }

  //
  const Pdb::Model::coords Pdb::Model::getMinCoords() const {
    float minX = std::numeric_limits<float>::max();
    float minY = std::numeric_limits<float>::max();
    float minZ = std::numeric_limits<float>::max();

    for (auto&& [x, y, z] : coords_) {
      if (x < minX) minX = x;
      if (y < minY) minY = y;
      if (z < minZ) minZ = z;
    }

    return std::make_tuple(minX, minY, minZ);
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
            fileCorrupted(path);
          }
          else {
            containsModels = true;
            awaitsENDMDL = true;
            std::string modelID = getKeyword(line, modelSerialFrom_, modelSerialTo_);
            if (modelID.empty()) {
              fileCorrupted(path);
            }
            models_.emplace_back(modelID);
          }
        }

        else if (section == "ENDMDL") {
          if (awaitsENDMDL) {
            awaitsENDMDL = false;
            ++currentModelIndex;
          }
          else {
            fileCorrupted(path);
          }
        }

        else if (section == "ATOM") {
          if (!awaitsHetTER) {
            if (!containsModels) {
              containsModels = true;
              models_.emplace_back("1");
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
              fileCorrupted(path);
            }
          }
          else {
            fileCorrupted(path);
          }
        }

        else if (section == "HETATM") {
          if (!awaitsAtomTER) {
            if (!containsModels) {
              containsModels = true;
              models_.emplace_back("1");
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
              fileCorrupted(path);
            }
          }
          else {
            fileCorrupted(path);
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
      fileCorrupted(path);
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
      if (!line[i]) break;
      else if (!isspace(line[i])) {
        keyword += line[i];
      }
    }

    return keyword;
  }

}

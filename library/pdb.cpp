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

  /*
    Fills a provided vector with Atom class instances belonging to a residue.
    Parameters:
    * atomVector: vector of Atom class instances, which will be filled with atoms belonging to a Residue class instance.
  */
  void Pdb::Model::Chain::Residue::getAtoms(Pdb::Model::Chain::Residue::atoms& atomVector) const {
    atomVector.clear();
    for (std::size_t i = indexFrom_; i <= indexTo_; ++i) {
      atomVector.emplace_back(model_, i);
    }
  }

  /*
    Finds and outputs an Atom class instance with a given ID.
    Parameters:
    * id: ID of the atom, must be present in Pdb::Model::atomSerials_ private variable.
    Returns: An Atom class instance representing the found atom.
  */
  const Pdb::Model::Chain::Residue::Atom Pdb::Model::Chain::Residue::findAtom(const std::string& id) {
    for (std::size_t i = indexFrom_; i <= indexTo_; ++i) {
      if (model_.atomSerials_[i] == id) {
        return Pdb::Model::Chain::Residue::Atom(model_, i);
      }
    }
    throw std::invalid_argument("Atom ID " + id + " was not found!");
  }

  /*
    Fills a provided vector with HetAtom class instances belonging to a ligand.
    Parameters:
    * atomVector: vector of HetAtom class instances, which will be filled with atoms belonging to a HetResidue class instance.
  */
  void Pdb::Model::Chain::HetResidue::getAtoms(Pdb::Model::Chain::HetResidue::hetAtoms& atomVector) const {
    atomVector.clear();
    for (std::size_t i = indexFrom_; i <= indexTo_; ++i) {
      atomVector.emplace_back(model_, i);
    }
  }

  /*
    Finds and outputs a HetAtom class instance with a given ID.
    Parameters:
    * id: ID of the atom, must be present in Pdb::Model::hetAtomSerials_ private variable.
    Returns: A HetAtom class instance representing the found atom.
  */
  const Pdb::Model::Chain::HetResidue::HetAtom Pdb::Model::Chain::HetResidue::findAtom(const std::string& id) {
    for (std::size_t i = indexFrom_; i <= indexTo_; ++i) {
      if (model_.hetAtomSerials_[i] == id) {
        return Pdb::Model::Chain::HetResidue::HetAtom(model_, i);
      }
    }
    throw std::invalid_argument("HET Atom ID " + id + " was not found!");
  }

  /*
    Counts all residues present within a chain.
    Returns: Number of residues present within a chain.
  */
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

  /*
    Counts all ligands present within a chain.
    Returns: Number of ligands present within a chain.
  */
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

  /*
    Counts all atoms belonging to ligands present within a chain.
    Returns: Number of atoms belonging to ligands present within a chain.
  */
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

  /*
    Fills a provided vector with Residue class instances belonging to a chain.
    Parameters:
    * residueVector: vector of Residue class instances, which will be filled with residues belonging to a Chain class instance.
  */
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

  /*
    Fills a provided vector with HetResidue class instances belonging to a chain.
    Parameters:
    * residueVector: vector of HetResidue class instances, which will be filled with ligands belonging to a Chain class instance.
  */
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

  /*
    Finds and outputs a Residue class instance with a given ID.
    Parameters:
    * id: ID of the residue, must be present in Pdb::Model::resSeqs_ private variable.
    Returns: A Residue class instance representing the found residue.
  */
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

  /*
    Finds and outputs a HetResidue class instance with a given ID.
    Parameters:
    * id: ID of the ligand, must be present in Pdb::Model::hetResSeqs_ private variable.
    Returns: A HetResidue class instance representing the found ligand.
  */
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

  /*
    Outputs the size of Pdb::Model::surfaceResidues_ private variable.
    If is has no content, a method is called to set value to it.
    Returns: Number of residues accessible to the solvent and number of residues inaccessible to the solvent.
  */
  const Pdb::Model::size_tTuple Pdb::Model::getNumberOfSurfaceAndBuried() {
    getSurfaceResidues();

    return std::make_tuple(
      surfaceResidues_.size(),
      getNumberOfResidues() - surfaceResidues_.size()
    );
  }

  /*
    Outputs how many times is each residue type present as solvent accessible and solvent inaccessible, by analyzing Pdb::Model::surfaceResidues_ private variable.
    If is has no content, a method is called to set value to it.
    Returns: Statistics of residue type count on the surface and in the core.
  */
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

  /*
    Outputs the percentage of polar residues in those accessible by solvent and buried in the core by analyzing Pdb::Model::surfaceStats_, Pdb::Model::buriedStats_ and Pdb::Model::surfaceResidues_.
    Returns: Portions of polar residues.
  */
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

  /*
    Outputs the distance between two furthest atoms in a model, from Pdb::Model::farthestAtoms_.
    If is has no content, a method is called to set value to it.
    Returns: Width of the model.
  */
  const float Pdb::Model::getWidth() {
    if (std::get<0>(farthestAtoms_) == 0) {
      getFarthestAtoms();
    }

    return std::get<0>(farthestAtoms_);
  }

  /*
    Outputs the diameter of the smallest circumsphere of all atoms in a model.
    If two farthest atoms have not been found yet, a method is called to find them.
    Then we set the mid location between them as the center of the circumsphere and try to find an atom outside the sphere.
    If such atom is not found, we return the width of the model.
    Else we find a new center as the circumcenter of the three atoms (two farthest and the one farthest from the old center) and make the radius slightly larger.
    Iterate until no atom outside the new sphere is found.
    Notice: This method return an approximate value!
    Returns: Diameter of models circumsphere.
  */
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

  /*
    Fills a provided vector with Atom class instances that are within a maximum distance from any atom of given ligand.
    Parameters:
    * het: a ligand, atoms within a maximum distance from its atoms shall be found.
    * atoms: vector of Atom class instances, which will be filled with atoms that are within a maximum distance from any atom of het.
    * maxDistance: maximum distance within which atoms will be searched.
  */
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

  /*
    Fills a provided vector with Residue class instances, at least one of its atoms have to be within a maximum distance from any atom of given ligand.
    Parameters:
    * het: a ligand, residues within a maximum distance from its atoms shall be found.
    * residues: vector of Residue class instances, which will be filled with residues that have at least one atom within a maximum distance from any atom of het.
    * maxDistance: maximum distance within which atoms will be searched.
  */
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

  /*
    Fills a provided vector with Chain class instances belonging to a model.
    Parameters:
    * chainVector: vector of Chain class instances, which will be filled with chains belonging to a Model class instance.
  */
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

  /*
    Finds and outputs a Chain class instance with a given ID.
    Parameters:
    * id: ID of the chain, must be present in Pdb::Model::chainIDs_ private variable.
    Returns: A Chain class instance representing the found chain.
  */
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

  /*
    Computes the circumcenter of a triangle with 3D coordinates.
    Parameters:
    * a: First point of the triangle.
    * b: Second point of the triangle.
    * c: Third point of the triangle.
    Returns: Coordinates of the computed circumcenter.
  */
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

  /*
    Sets Pdb::Model::farthestAtoms_ to contain the two farthest atoms in a model and their distance, if the values have not been set previously.
    It first creates a convex hull of the atoms, if it hasn't been created previously.
    Then it applies the naive algorithm to find the most distant pair of atoms.
    Time complexity is O(n log n + h^2) in the worst case, where n is the number of atoms in the model and h is the number of atoms in the convex hull.
    The convex hull should contain only a fraction of all atoms, therefore being more efficient than applying the naive O(n^2) algorithm to the whole model.
  */
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

  /*
    Uses the QuickHull library to compute the convex hull of a model and saves it in Pdb::Model::convexHull_ private variable, if it hasn't been computed previously.
    Time complexity of the quickhull algorithm is O(n log n), n is the number of atoms.
    Returns: True iff this method is called for the first time.
  */
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

  /*
    Checks if the name of given atom does not equal any name associated with the backbone.
    Parameters:
    * atom: An instance of Atom class.
    Returns: True iff the atom name is not associated with the backbone.
  */
  const bool Pdb::Model::atomIsSideChain(const Chain::Residue::Atom& atom) {
    const std::string& atomName = atom.getAtomName();

    if (atomName == "N") return false;
    if (atomName == "CA") return false;
    if (atomName == "C") return false;
    if (atomName == "O") return false;
    if (atomName == "OXT") return false;
    return true;
  }

  /*
    Stores residues with at least one atom accessible to a solvent to Pdb::Model::surfaceResidues_ private variable.
    First we build a 3D grid having size slighly larger than the whole model.
    One cell of the grid is just as large as the solvent molecule.
    Then we place all atoms of the model into the grid based on the following:
    If the atoms has coordinates (x, y, z), then the atom will be placed into the grid on position [x][y][z].
    Obviously, the coordinates have to be transformed, as arrays do not accept negative indexes and the grid is scaled to the solvent size.
    Notice: More atoms can belong into one cell of the grid, but only one can be stored.
    Notice: Such atoms belonging into one cell are usually part of the same residue.
    At the end, we traverse empty cells of the 3D grid with BFS as if we were scanning the surroundings of the model with a solvent molecule.
    This yields approximately the solvent accessible residues.
  */
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

  /*
    Finds the minimum value of each coordinate among atoms in a model.
    This method is used to fetch the offset from [0][0][0] cell of the 3D grid from getSurfaceResidues().
    Returns: Triple of minimum coordinate values.
  */
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

  /*
    Finds and outputs a Model class instance with a given ID.
    Parameters:
    * id: ID of the model, must be present in Pdb::models_ private variable.
    Returns: A Model class instance representing the found residue.
  */
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

  /*
    Goes through each line of a .pdb file and enters parsed data into Pdb abd Model class intances' private variables.
    Parameters:
    * path: path to a .pdb file which shall be parsed
  */
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

  /*
    Outputs a part of given line from a given position to a given position.
    Parameters:
    * line: one line of a .pdb file.
    * from: position on line from which the part of line will be parsed.
    * to: position on line to which the part of line will be parsed.
    Returns: A part of line.
  */
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

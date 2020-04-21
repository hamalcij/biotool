// biotool
// Bioinformatics Toolbox
//
// library/pdb.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// Definition of class Pdb
//

#ifndef PDB_HPP
#define PDB_HPP

#include<string>
#include<vector>
#include<set>
#include<tuple>
#include<fstream>
#include<exception>
#include<cctype>

////////////////////////////////////////////////////////////////////////
// BEGIN PDB

class Pdb {
public:

  ////////////////////////////////////////////////////////////////////////
  // BEGIN MODEL

  class Model {
  public:

    using coords = std::tuple<float, float, float>;
    using coordinates = std::vector<coords>;
    using chars = std::vector<char>;
    using floats = std::vector<float>;
    using ints = std::vector<std::string>;
    using strings = std::vector<std::string>;

    ////////////////////////////////////////////////////////////////////////
    // BEGIN CHAIN

    class Chain {
    public:

      ////////////////////////////////////////////////////////////////////////
      // BEGIN RESIDUE

      class Residue {
      public:

        ////////////////////////////////////////////////////////////////////////
        // BEGIN ATOM

        class Atom {
        public:
          Atom(
            const Model& model,
            const std::size_t index
          ) :
            model_(model),
            index_(index)
          {}

          const std::string& getID() const { return model_.atomSerials_[index_]; }
          const std::string& getAtomName() const { return model_.atomNames_[index_]; }
          const coords& getCoords() const { return model_.coords_[index_]; }
          const float& getOcuupancy() const { return model_.occupancies_[index_]; }
          const float& getTempFactor() const { return model_.tempFactors_[index_]; }
          const std::string& getElement() const { return model_.elements_[index_]; }
          const std::string& getCharge() const { return model_.charges_[index_]; }

        private:
          const Model& model_;
          const std::size_t index_;
        };

        // END ATOM
        ////////////////////////////////////////////////////////////////////////

        using atoms = std::vector<Atom>;

        Residue(
          const Model& model,
          const std::size_t indexFrom,
          const std::size_t indexTo
        ) :
          model_(model),
          indexFrom_(indexFrom),
          indexTo_(indexTo)
        {}

        const std::string& getID() const { return model_.resSeqs_[indexFrom_]; }
        const std::string& getResidueName() const { return model_.resNames_[indexFrom_]; }
        const std::size_t getNumberOfAtoms() const { return indexTo_ - indexFrom_ + 1; }

        void getAtoms(atoms& atomVector) const;
        const Atom findAtom(const std::string& id);

      private:
        const Model& model_;
        const std::size_t indexFrom_;
        const std::size_t indexTo_;
      };

      // END RESIDUE
      ////////////////////////////////////////////////////////////////////////
      // BEGIN HETRESIDUE

      class HetResidue {
      public:

        ////////////////////////////////////////////////////////////////////////
        // BEGIN HETATOM

        class HetAtom {
        public:
          HetAtom(
            const Model& model,
            const std::size_t index
          ) :
            model_(model),
            index_(index)
          {}

          const std::string& getID() const { return model_.hetAtomSerials_[index_]; }
          const std::string& getAtomName() const { return model_.hetAtomNames_[index_]; }
          const coords& getCoords() const { return model_.hetCoords_[index_]; }
          const float& getOcuupancy() const { return model_.hetOccupancies_[index_]; }
          const float& getTempFactor() const { return model_.hetTempFactors_[index_]; }
          const std::string& getElement() const { return model_.hetElements_[index_]; }
          const std::string& getCharge() const { return model_.hetCharges_[index_]; }

        private:
          const Model& model_;
          const std::size_t index_;
        };

        // END HETATOM
        ////////////////////////////////////////////////////////////////////////

        using hetAtoms = std::vector<HetAtom>;

        HetResidue(
          const Model& model,
          const std::size_t indexFrom,
          const std::size_t indexTo
        ) :
          model_(model),
          indexFrom_(indexFrom),
          indexTo_(indexTo)
        {}

        const std::string& getID() const { return model_.hetResSeqs_[indexFrom_]; }
        const std::string& getResidueName() const { return model_.hetResNames_[indexFrom_]; }
        const std::size_t getNumberOfHetAtoms() const { return indexTo_ - indexFrom_ + 1; }

        void getAtoms(hetAtoms& atomVector) const;
        const HetAtom findAtom(const std::string& id);

      private:
        const Model& model_;
        const std::size_t indexFrom_;
        const std::size_t indexTo_;
      };

      // END HETRESIDUE
      ////////////////////////////////////////////////////////////////////////

      using residues = std::vector<Residue>;
      using hetResidues = std::vector<HetResidue>;

      Chain(
        const Model& model,
        const std::size_t indexFrom,
        const std::size_t indexTo
      ) :
        model_(model),
        indexFrom_(indexFrom),
        indexTo_(indexTo)
      {}

      const char getID() const { return model_.chainIDs_[indexFrom_]; }
      const std::size_t getNumberOfResidues() const;
      const std::size_t getNumberOfAtoms() const { return indexTo_ - indexFrom_ + 1; }
      const std::size_t getNumberOfHetResidues() const;
      const std::size_t getNumberOfHetAtoms() const;

      void getResidues(residues& residueVector) const;
      void getHetResidues(hetResidues& residueVector) const;

      const Residue findResidue(const std::string& id) const;
      const HetResidue findHetResidue(const std::string& id) const;

    private:
      const Model& model_;
      const std::size_t indexFrom_;
      const std::size_t indexTo_;
    };

    // END CHAIN
    ////////////////////////////////////////////////////////////////////////

    using chains = std::vector<Chain>;

    friend class Pdb;

    Model(Pdb& pdb, const std::string& id) : pdb_(pdb), id_(id) {}

    const std::string& getID() const { return id_; }
    const std::size_t getNumberOfChains() const {
      return getNumberOfDistinctElements(chainIDs_);
    }
    const std::size_t getNumberOfResidues() const {
      return getNumberOfDistinctElements(resSeqs_);
    }
    const std::size_t getNumberOfHetResidues() const {
      return getNumberOfDistinctElements(hetResSeqs_);
    }
    const std::size_t getNumberOfAtoms() const { return atomSerials_.size(); }
    const std::size_t getNumberOfHetAtoms() const { return hetAtomSerials_.size(); }
    const std::size_t getWidth() const;
    const Chain::Residue::atoms getAtomsCloseToLigand() const;
    const Chain::residues getResiduesCloseToLigand() const;

    void getChains(chains& chainVector) const;
    const Chain findChain(const char id) const;

  private:
    template<typename T>
    const std::size_t getNumberOfDistinctElements(const std::vector<T>& column) const {
      if (column.empty()) {
        return 0;
      }
      std::size_t numberOfDistinctElements = 1;
      T current = column.front();
      for (auto&& data : column) {
        if (current != data) {
          current = data;
          ++numberOfDistinctElements;
        }
      }
      return numberOfDistinctElements;
    }

    Pdb& pdb_;
    const std::string id_;

    ints atomSerials_;
    strings atomNames_;
    chars altLocs_;
    strings resNames_;
    chars chainIDs_;
    ints resSeqs_;
    chars iCodes_;
    coordinates coords_;
    floats occupancies_;
    floats tempFactors_;
    strings elements_;
    strings charges_;

    ints hetAtomSerials_;
    strings hetAtomNames_;
    chars hetAltLocs_;
    strings hetResNames_;
    chars hetChainIDs_;
    ints hetResSeqs_;
    chars hetICodes_;
    coordinates hetCoords_;
    floats hetOccupancies_;
    floats hetTempFactors_;
    strings hetElements_;
    strings hetCharges_;
  };

  // END MODEL
  ////////////////////////////////////////////////////////////////////////

  using models = std::vector<Model>;

  Pdb(const std::string& path) { parsePdbFile(path); }

  const std::size_t getNumberOfModels() const { return models_.size(); }
  const models& getModels() const { return models_; }
  const Model& findModel(const std::string& id) const;

private:
  void parsePdbFile(const std::string& path);
  std::string getKeyword(char* line, const unsigned short from, const unsigned short to);
  void pdbFileCorrupted(const std::string& path);
  bool fileIsEmpty(std::ifstream& file);

  models models_;

  static const unsigned short recordNameFrom_{0};
  static const unsigned short recordNameTo_{5};
  static const unsigned short modelSerialFrom_{10};
  static const unsigned short modelSerialTo_{13};
  static const unsigned short atomSerialFrom_{6};
  static const unsigned short atomSerialTo_{10};
  static const unsigned short atomNameFrom_{12};
  static const unsigned short atomNameTo_{15};
  static const unsigned short altLocPos_{16};
  static const unsigned short resNameFrom_{17};
  static const unsigned short resNameTo_{19};
  static const unsigned short chainIDPos_{21};
  static const unsigned short resSeqFrom_{22};
  static const unsigned short resSeqTo_{25};
  static const unsigned short iCodePos_{26};
  static const unsigned short xFrom_{30};
  static const unsigned short xTo_{37};
  static const unsigned short yFrom_{38};
  static const unsigned short yTo_{45};
  static const unsigned short zFrom_{46};
  static const unsigned short zTo_{53};
  static const unsigned short occupancyFrom_{54};
  static const unsigned short occupancyTo_{59};
  static const unsigned short tempFactorFrom_{60};
  static const unsigned short tempFactorTo_{65};
  static const unsigned short elementFrom_{76};
  static const unsigned short elementTo_{77};
  static const unsigned short chargeFrom_{78};
  static const unsigned short chargeTo_{79};
  static const unsigned short bufferSize_{82};
};

// END PDB
////////////////////////////////////////////////////////////////////////

#endif // PDB_HPP
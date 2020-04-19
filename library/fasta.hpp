// biotool
// Bioinformatics Toolbox
//
// library/fasta.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// Definition of class Fasta
//

#ifndef FASTA_HPP
#define FASTA_HPP

#include<vector>
#include<string>
#include<fstream>
#include<exception>

class Fasta {
public:

  using descriptions = std::vector<std::string>;
  using sequences = std::vector<std::string>;

  class Molecule {
  public:
    Molecule(Fasta& fasta, const std::size_t& index) : fasta_(fasta), index_(index) {}

    std::string getDescription() { return fasta_.descriptions_[index_]; }
    std::string getSequence() { return fasta_.sequences_[index_]; }
    std::string getSubSequence(const std::size_t& from, const std::size_t& to) { return fasta_.sequences_[index_].substr(from - 1, to - from + 1); }
    std::size_t getSequenceLength() { return fasta_.sequences_[index_].length(); }

  private:
    Fasta& fasta_;
    const std::size_t index_;
  };

  Fasta(const std::string& path) { parseFastaFile(path); }

  const std::size_t numberOfMolecules() const { return descriptions_.size(); }

private:
  void parseFastaFile(const std::string& path);
  void fastaFileCorrupted(const std::string& path);
  bool fileIsEmpty(std::ifstream& file);

  descriptions descriptions_;
  sequences sequences_;
};

#endif // FASTA_HPP

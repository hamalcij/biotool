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
// library/fasta.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// Definition of class Fasta
//

#ifndef FASTA_HPP
#define FASTA_HPP

#include "file.hpp"

#include<vector>

namespace biotool {

  /*
    Parses a .fasta file
  */
  class Fasta {
  public:

    using descriptions = std::vector<std::string>;
    using sequences = std::vector<std::string>;

    /*
      Proxy class to access parsed data in an instance of Fasta class
    */
    class Molecule {
    public:

      /*
        Molecule class ctor
        Parameters:
        * fasta: instance of Fasta class, will be accessed
        * index: index in the database of parsed .fasta file in fasta
      */
      Molecule(Fasta& fasta, const std::size_t& index) : fasta_(fasta), index_(index) {}

      std::string getDescription() { return fasta_.descriptions_[index_]; }
      std::string getSequence() { return fasta_.sequences_[index_]; }
      std::string getSubSequence(const std::size_t& from, const std::size_t& to) { return fasta_.sequences_[index_].substr(from - 1, to - from + 1); }
      std::size_t getSequenceLength() { return fasta_.sequences_[index_].length(); }

    private:
      Fasta& fasta_;
      const std::size_t index_;
    };

    /*
      Fasta class ctor
      Parameters:
      * path: path to a .fasta file which shall be parsed
    */
    Fasta(const std::string& path) { parseFastaFile(path); }

    const std::size_t numberOfMolecules() const { return descriptions_.size(); }

  private:
    void parseFastaFile(const std::string& path);

    descriptions descriptions_;
    sequences sequences_;
  };

}

#endif // FASTA_HPP

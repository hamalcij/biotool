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
// library/sequencePair.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// Definition of Hamming class
//

#ifndef SEQUENCE_PAIR_HPP
#define SEQUENCE_PAIR_HPP

#include<string>
#include<vector>
#include<exception>
#include<tuple>

namespace biotool {

  /*
    Stores two sequences and provides methods to analyse them
  */
  class SequencePair {
  public:

    using sequences = std::tuple<std::string, std::string>;
    using optimalAlignments = std::vector<sequences>;
    using alignment = std::tuple<std::size_t, optimalAlignments>;

    /*
      SequencePair class ctor
      Parameters:
      * seq1: one sequence, usually of a protein or nucleic acid
      * seq2: another sequence, usually of a protein or nucleic acid
    */
    SequencePair(const std::string& seq1, const std::string& seq2) : sequence1_(seq1), sequence2_(seq2) {}

    const sequences getSequences() const { return std::make_tuple(sequence1_, sequence2_); }

    const std::size_t getHammingDistance() const;
    const alignment getEditDistance() const;

  private:
    void checkSequenceLength() const;
    const std::size_t dpMin(const std::size_t& insertion, const std::size_t& deletion, const std::size_t& substitution) const;

    const std::string sequence1_;
    const std::string sequence2_;
  };

}

#endif // SEQUENCE_PAIR_HPP

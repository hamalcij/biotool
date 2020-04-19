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

using sequences = std::tuple<std::string, std::string>;

class SequencePair {
public:
  SequencePair(const std::string& seq1, const std::string& seq2) : sequence1_(seq1), sequence2_(seq2) {}

  const sequences getSequences() const { return std::make_tuple(sequence1_, sequence2_); }

  const std::size_t getHammingDistance() const;

private:
  void checkSequenceLength() const;

  const std::string sequence1_;
  const std::string sequence2_;
};

#endif // SEQUENCE_PAIR_HPP

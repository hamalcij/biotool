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

class SequencePair {
public:

  using sequences = std::tuple<std::string, std::string>;
  using optimalAlignments = std::vector<sequences>;
  using alignment = std::tuple<std::size_t, optimalAlignments>;

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

#endif // SEQUENCE_PAIR_HPP

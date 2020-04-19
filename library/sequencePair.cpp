// biotool
// Bioinformatics Toolbox
//
// library/sequencePair.cpp
// Copyright (c) 2020 Hamalčík Jan
//
// Implementation of Hamming class
//

#include "sequencePair.hpp"

const std::size_t SequencePair::getHammingDistance() const {
  checkSequenceLength();

  std::size_t hammingDistance = 0;
  for (std::size_t i = 0; i < sequence1_.length(); ++i) {
    if (sequence1_[i] != sequence2_[i]) {
      ++hammingDistance;
    }
  }

  return hammingDistance;
}

void SequencePair::checkSequenceLength() const {
  if (sequence1_.length() != sequence2_.length()) {
    throw std::length_error("Hamming distance can only be calculated for sequences with same length!");
  }
}

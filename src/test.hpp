// biotool
// Bioinformatics Toolbox
//
// src/test.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// Test showing how biotool API works
//

#include "../include/biotool.hpp"

#include<iostream>

class TestFasta {
public:
  TestFasta(const std::string& file) : myFasta_(FASTA(file)) {}

  void performTest() {
    for (std::size_t i = 0; i < NUMBER_OF_MOLECULES(myFasta_); ++i) {
      MOLECULE myMolecule(myFasta_, i);
      std::cout << "NEW MOLECULE:" << std::endl;
      std::cout << "Description: " << DESCRIPTION(myMolecule) << std::endl;
      std::cout << "Sequence: " << SEQUENCE(myMolecule) << std::endl;
      std::cout << "Sequence length: "<< SEQUENCE_LENGTH(myMolecule) << std::endl;;
      std::cout << "Subsequence 2...24: " << SUBSEQUENCE(myMolecule, 2, 24) << std::endl;
      std::cout << std::endl;
    }
  }

private:
  FASTA myFasta_;
};

class TestSequencePair {
public:
  TestSequencePair(const std::string& seq1, const std::string& seq2) : mySequencePair_(SEQUENCE_PAIR(seq1, seq2)) {}

  void performTest() {
    auto[sequence1, sequence2] = SEQUENCES(mySequencePair_);

    std::cout << "Sequence 1: " << sequence1 << std::endl;
    std::cout << "Sequence 2: " << sequence2 << std::endl;
    std::cout << "Hamming distance: " << HAMMING_DISTANCE(mySequencePair_) << std::endl;
  }

private:
  SEQUENCE_PAIR mySequencePair_;
};

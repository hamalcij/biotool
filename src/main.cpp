// biotool
// Bioinformatics Toolbox
//
// src/main.cpp
// Copyright (c) 2020 Hamalčík Jan
//
// Test showing how biotool API works
//

#include "test.hpp"

int main() {
  TestFasta fasta("/Users/janhamalcik/github/biotool/testfiles/fasta");
  fasta.performTest();

  TestSequencePair sequencePair("SHORTS", "HORSTS");
  sequencePair.performTest();

  return 0;
}

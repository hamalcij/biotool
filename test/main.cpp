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
// test/main.cpp
// Copyright (c) 2020 Hamalčík Jan
//
// Test showing how biotool API works
//

#include "test.hpp"

int main() {
  std::cout << "This is biotool's test program." << std::endl;

  for (;;) {
    int answer;
    std::string path;
    std::string seq1;
    std::string seq2;

    std::cout << std::endl;
    std::cout << "########################################################" << std::endl;
    std::cout << std::endl;
    std::cout << "To test the FASTA parser, type: 1" << std::endl;
    std::cout << "To test the sequence pair analysis, type: 2" << std::endl;
    std::cout << "To test the PDB parser, type: 3" << std::endl;
    std::cout << "To test the Clustal parser, type: 4" << std::endl;
    std::cout << "To exit biotool's test program, type: 0" << std::endl;
    std::cin >> answer;
    std::cout << std::endl;

    if (answer == 1) {
      std::cout << "Type full path to your .fasta file:" << std::endl;
      std::cin >> path;
      TestFasta fasta(path);
      fasta.performTest();
    }

    else if (answer == 2) {
      std::cout << "Type first sequence:" << std::endl;
      std::cin >> seq1;
      std::cout << "Type second sequence:" << std::endl;
      std::cin >> seq2;
      TestSequencePair sequencePair(seq1, seq2);
      sequencePair.performTest();
    }

    else if (answer == 3) {
      std::cout << "Type full path to your .pdb file:" << std::endl;
      std::cin >> path;
      TestPdb pdb(path);
      pdb.performTest();
    }

    else if (answer == 4) {
      std::cout << "Type full path to your clustal file:" << std::endl;
      std::cin >> path;
      TestClustal clustal(path);
      clustal.performTest();
    }

    else if (answer == 0) {
      return 0;
    }

    else {
      throw std::invalid_argument("\"" + std::to_string(answer) + "\" is an invalid argument!");
    }
  }
}

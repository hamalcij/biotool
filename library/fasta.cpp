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
// library/fasta.cpp
// Copyright (c) 2020 Hamalčík Jan
//
// Implementation of class Fasta
//

#include "fasta.hpp"

namespace biotool {

  /*
    Goes through each line of a .fasta file and enters parsed data to descriptions_ and sequences_ private variables.
    Parameters:
    * path: path to a .fasta file which shall be parsed
  */
  void Fasta::parseFastaFile(const std::string& path) {
    std::ifstream file(path);

    if (file.good() && !fileIsEmpty(file)) {
      std::string line;
      std::string sequence = "";
      bool sequenceIsBeingRead = false;
      bool sequenceWasRead = true;

      while(std::getline(file, line)) {
        if (line[0] == '>') {
          if (sequenceWasRead) {
            sequenceIsBeingRead = false;
            sequenceWasRead = false;
            descriptions_.emplace_back(line.substr(1));
          }
          else {
            fileCorrupted(path);
          }
        }
        else if (line.empty()) {
          if (!sequenceWasRead) {
            fileCorrupted(path);
          }
        }
        else {
          if (sequenceIsBeingRead) {
            sequences_.back() += line;
          }
          else {
            sequenceIsBeingRead = true;
            sequenceWasRead = true;
            sequences_.push_back(line);
          }
        }
      }
    }
    else {
      fileCorrupted(path);
    }

    file.close();
  }

}

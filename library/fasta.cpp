// biotool
// Bioinformatics Toolbox
//
// library/fastaParser.cpp
// Copyright (c) 2020 Hamalčík Jan
//
// Implementation of class Fasta
//

#include "fasta.hpp"

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
          fastaFileCorrupted(path);
        }
      }
      else if (line.empty()) {
        if (!sequenceWasRead) {
          fastaFileCorrupted(path);
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
    fastaFileCorrupted(path);
  }

  file.close();
}

void Fasta::fastaFileCorrupted(const std::string& path) {
  throw std::ifstream::failure("Fasta file at " + path + " is corrupted!");

}

bool Fasta::fileIsEmpty(std::ifstream& file) {
  return file.peek() == std::ifstream::traits_type::eof();
}

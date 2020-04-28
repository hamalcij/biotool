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
// library/clustal.hpp
// Copyright (c) 2020 Hamalčík Jan
//
// Definition of class Clustal
//

#include "clustal.hpp"

namespace biotool {

  /*
    Outputs a sequence of MSA with described with given ID.
    * id: Entry ID of MSA.
    Returns: Sequence described with id.
  */
  const std::string Clustal::getSequence(const std::string& id) const {
    auto it = std::find(ids_.cbegin(), ids_.cend(), id);
    if (it == ids_.cend()) {
      throw std::invalid_argument("Sequence ID " + id + " was not found!");
    }
    std::size_t i = it - ids_.cbegin();

    return readSequenceByIndex(i);
  }

  /*
    Goes through each line of a clustal file and enters parsed data into Clustal class instance's private variables.
    Parameters:
    * path: path to a .pdb file which shall be parsed
  */
  void Clustal::parseClustalFile(const std::string& path) {
    std::ifstream file(path);

    if (file.good() && !fileIsEmpty(file)) {
      char line[bufferSize_];
      bool isInsindeSequenceBlock = false;
      bool isFirstBlock = true;
      std::size_t i = 0;
      std::size_t sequenceStartIndex = 0;
      std::size_t sequenceEndIndex = 0;
      std::size_t emplaceFromIndex = 0;
      std::size_t columnIndex = 0;
      std::string id = "";

      file.getline(line, bufferSize_);
      id = readAWord(line, i);
      if (file.fail() || (id != "CLUSTAL" && id != "CLUSTALW")) {
        fileCorrupted(path);
      }

      try {
        for (;;) {
          file.getline(line, bufferSize_);
          if (file.fail()) break;

          i = 0;
          if (line[i]) {

            if (isspace(line[i])) {

              // Conservation line
              if (isInsindeSequenceBlock) {
                i = sequenceStartIndex;
                while (line[i] && i <= sequenceEndIndex) {
                  conservation_.emplace_back(line[i]);
                  ++i;
                }
                if (i < sequenceEndIndex) {
                  while (i <= sequenceEndIndex) {
                    conservation_.emplace_back(' ');
                    ++i;
                  }
                }
                isInsindeSequenceBlock = false;
                isFirstBlock = false;
              }
            }

            else {
              if (!isInsindeSequenceBlock) {
                emplaceFromIndex = alignment_.size();
              }
              id = readAWord(line, i);
              while (isspace(line[i])) {
                ++i;
              }
              sequenceStartIndex = i;
              if (!isInsindeSequenceBlock) {
                columnIndex = 0;
                if (isFirstBlock) {
                  ids_.emplace_back(id);
                }
                else {
                  if (ids_[columnIndex] != id) throw std::invalid_argument("");
                }
                while (line[i] && !isspace(line[i])) {
                  alignment_.emplace_back(
                    std::make_shared<Clustal::column>(1, toupper(line[i]))
                  );
                  ++i;
                }
                sequenceEndIndex = i - 1;
                isInsindeSequenceBlock = true;
              }
              else {
                ++columnIndex;
                if (!isFirstBlock) {
                  if (ids_[columnIndex] != id) throw std::invalid_argument("");
                }
                else {
                  ids_.emplace_back(id);
                }
                while (line[i] && !isspace(line[i])) {
                  alignment_[emplaceFromIndex + i - sequenceStartIndex]
                  ->emplace_back(toupper(line[i]));
                  ++i;
                }
              }
            }
          }
          else {
            if (isInsindeSequenceBlock) {
              for (std::size_t j = sequenceStartIndex; j <= sequenceEndIndex; ++j) {
                conservation_.emplace_back(' ');
              }
              isInsindeSequenceBlock = false;
              isFirstBlock = false;
            }
          }
        }
      }
      catch (...) {
        fileCorrupted(path);
      }
    }
    else {
      fileCorrupted(path);
    }

    file.close();
  }

  /*
    Outputs a part of given line from given index until the first space or end of line.
    Parameters:
    * line: one line of a clustal file.
    * i: position on line from which the part of line will be parsed.
    Returns: A part of line.
  */
  const std::string Clustal::readAWord(const char* line, std::size_t& i) {
    std::string word = "";
    while (line[i] && !isspace(line[i])) {
      word += line[i];
      ++i;
    }

    return word;
  }

  /*
    Outputs a sequence of MSA on given line within a block.
    * i: index of a line within a block.
    Returns: Sequence on i'th line within a block.
  */
  const std::string Clustal::readSequenceByIndex(const std::size_t i) const {
    std::string sequence = "";
    for (auto&& column : alignment_) {
      sequence += (*column)[i];
    }

    return sequence;
  }

}

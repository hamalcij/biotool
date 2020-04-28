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
// library/sequencePair.cpp
// Copyright (c) 2020 Hamalčík Jan
//
// Implementation of Hamming class
//

#include "sequencePair.hpp"

namespace biotool {

  /*
    Computes Hamming distance of private variables sequence1_ and sequence2_, if they have the same length.
    Returns: The Hamming distance.
  */
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

  /*
    Computes edit distance of private variables sequence1_ and sequence2_ and recreates all optimal alignments.
    Dynamic programming is used, each penalty is scored 1, backtracking is used to obtain all optimal alignments.
    Returns: The edit distance and all alignments yielding the same edit distance.
  */
  const SequencePair::alignment SequencePair::getEditDistance() const {
    const std::size_t sequence1Size = sequence1_.length() + 1;
    const std::size_t sequence2Size = sequence2_.length() + 1;

    std::size_t dpMatrix[sequence1Size][sequence2Size];

    for (std::size_t i = 0; i < sequence1Size; ++i) {
      dpMatrix[i][0] = i;
    }
    for (std::size_t j = 1; j < sequence2Size; ++j) {
      dpMatrix[0][j] = j;
    }

    for (std::size_t i = 1; i < sequence1Size; ++i) {
      for (std::size_t j = 1; j < sequence2Size; ++j) {
        if (sequence1_[i - 1] == sequence2_[j - 1]) {
          dpMatrix[i][j] = dpMatrix[i - 1][j - 1];
        }
        else {
          dpMatrix[i][j] = dpMin(dpMatrix[i - 1][j] + 1, dpMatrix[i][j - 1] + 1, dpMatrix[i - 1][j - 1] + 1);
        }
      }
    }

    using backtrackedAlignment = std::tuple<std::size_t, std::size_t, std::string, std::string>;
    using backtrackedAlignments = std::vector<backtrackedAlignment>;

    #define I std::get<0>(*it)
    #define J std::get<1>(*it)
    #define SEQ1 std::get<2>(*it)
    #define SEQ2 std::get<3>(*it)

    optimalAlignments alignments;
    backtrackedAlignments backtracking;
    backtracking.emplace_back(std::make_tuple(sequence1Size - 1, sequence2Size - 1, std::string(""), std::string("")));

    while (!backtracking.empty()) {
      backtrackedAlignments newBacktracking;

      for (auto it = backtracking.begin(); it != backtracking.end();) {
        if (I == 0) {
          while (J != 0) {
            --J;
            SEQ1 = "-" + SEQ1;
            SEQ2 = sequence2_[J] + SEQ2;
          }
          alignments.emplace_back(std::make_tuple(SEQ1, SEQ2));
          it = backtracking.erase(it);
        }
        else if (J == 0) {
          while (I != 0) {
            --I;
            SEQ1 = sequence1_[I] + SEQ1;
            SEQ2 = "-" + SEQ2;
          }
          alignments.emplace_back(std::make_tuple(SEQ1, SEQ2));
          it = backtracking.erase(it);
        }
        else {
          bool minNotIdentified = true;
          std::size_t min = dpMin(dpMatrix[I - 1][J], dpMatrix[I][J - 1], dpMatrix[I - 1][J - 1]);
          std::size_t thisI = I;
          std::size_t thisJ = J;
          std::string thisSeq1 = SEQ1;
          std::string thisSeq2 = SEQ2;

          if (min == dpMatrix[I - 1][J - 1]) {
            minNotIdentified = false;
            --I;
            --J;
            SEQ1 = sequence1_[I] + SEQ1;
            SEQ2 = sequence2_[J] + SEQ2;
          }
          if (min == dpMatrix[thisI - 1][thisJ]) {
            if (minNotIdentified) {
              minNotIdentified = false;
              --I;
              SEQ1 = sequence1_[I] + SEQ1;
              SEQ2 = "-" + SEQ2;
            }
            else {
              newBacktracking.emplace_back(std::make_tuple(thisI - 1, thisJ, sequence1_[thisI - 1] + thisSeq1, "-" + thisSeq2));
            }
          }
          if (min == dpMatrix[thisI][thisJ - 1]) {
            if (minNotIdentified) {
              --J;
              SEQ1 = "-" + SEQ1;
              SEQ2 = sequence2_[J] + SEQ2;
            }
            else {
              newBacktracking.emplace_back(std::make_tuple(thisI, thisJ - 1, "-" + thisSeq1, sequence2_[thisJ - 1] + thisSeq2));
            }
          }
          ++it;
        }
      }

      for (auto&& x : newBacktracking) {
        backtracking.push_back(x);
      }
    }

    return std::make_tuple(dpMatrix[sequence1Size - 1][sequence2Size - 1], alignments);
  }

  /*
    Checks if the private variables sequence1_ and sequence2_ have the same length.
    Used to check if Hamming distance may be computed.
  */
  void SequencePair::checkSequenceLength() const {
    if (sequence1_.length() != sequence2_.length()) {
      throw std::length_error("Hamming distance can only be calculated for sequences with same length!");
    }
  }

  /*
    Checks which of the 3 values considered in 2D dynamic programming and backtracking is the smallest.
    Parameters:
    * insertion: value of the cell at position [I-1][J] of the DP matrix.
    * deletion: value of the cell at position [I][J-1] of the DP matrix.
    * substitution: value of the cell at position [I-1][J-1] of the DP matrix.
    Returns: The lowest value of insertion, deletion and substitution.
  */
  const std::size_t SequencePair::dpMin(const std::size_t& insertion, const std::size_t& deletion, const std::size_t& substitution) const {
    std::size_t min = insertion;
    if (min > deletion) {
      min = deletion;
    }
    if (min > substitution) {
      min = substitution;
    }

    return min;
  }

}

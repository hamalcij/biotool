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

#ifndef CLUSTAL_HPP
#define CLUSTAL_HPP

#include "file.hpp"

#include<vector>
#include<memory>
#include<tuple>
#include<cctype>
#include<utility>
#include<algorithm>

namespace biotool {

  class Clustal {
  public:

    using column = std::vector<char>;
    using conservation = std::vector<char>;
    using ids = std::vector<std::string>;
    using alignment = std::vector<std::unique_ptr<column>>;

    Clustal(const std::string& path) { parseClustalFile(path); }

    const std::size_t getNumberOfSequences() const { return ids_.size(); }
    const std::size_t getNumberOfColumns() const { return alignment_.size(); }
    const std::string getSequence(const std::size_t i) const { return readSequenceByIndex(i); }
    const std::string getSequence(const std::string& id) const;
    const std::string& getID(const std::size_t i) const { return ids_[i]; }
    const column& getColumn(const std::size_t i) const { return *alignment_[i]; }

    template<typename Matrix, typename Order>
    const double sumOfPairs(const Matrix& matrix, const Order& order, const std::size_t index) const {
      double sum = 0;
      const column& col = *alignment_[index];

      for (std::size_t i = 0; i < col.size(); ++i) {
        for (std::size_t j = i + 1; j < col.size(); ++j) {
          if (isalpha(col[i]) && isalpha(col[j])) {
            const std::size_t iPos = findAminoAcidPosition(order, col[i]);
            const std::size_t jPos = findAminoAcidPosition(order, col[j]);
            if (matrix[iPos].size() > jPos) {
              sum += matrix[iPos][jPos];
            }
            else {
              sum += matrix[jPos][iPos];
            }
          }
        }
      }

      return sum;
    }

    template<typename Matrix, typename Order>
    const double sumOfPairs(const Matrix& matrix, const Order& order) const {
      double sum = 0;

      for (std::size_t i = 0; i < alignment_.size(); ++i) {
        sum += sumOfPairs(matrix, order, i);
      }

      return sum;
    }

    template<typename Order>
    const std::size_t findAminoAcidPosition(const Order& order, const char aminoAcid) const {
      for (std::size_t i = 0; i < order.size(); ++i) {
        if (order[i] == aminoAcid) return i;
      }
      throw std::domain_error("Amino acid " + std::string(1, aminoAcid) + " not found in scoring matrix!");
    }

  private:
    void parseClustalFile(const std::string& path);
    const std::string readAWord(const char* line, std::size_t& i);
    const std::string readSequenceByIndex(const std::size_t i) const;

    alignment alignment_;
    ids ids_;
    conservation conservation_;

    static const unsigned short bufferSize_{82};
  };

}

#endif // CLUSTAL_HPP

[//]: # (biotool)
[//]: # (Bioinformatics Toolbox)
[//]: # ()
[//]: # (docs/api/clustal.md)
[//]: # (Copyright © 2020 Hamalčík Jan)
[//]: # ()
[//]: # (Documentation of Clustal class API)
[//]: # ()

# Clustal API Documentation

Following macros will conveniently let you use `biotool`'s Clustal API:

## `CLUSTAL`

`CLUSTAL` macro defines a class that will let you access the Clustal
parser's methods, such as getting the best scoring columns stored in your
clustal file.

To initiate the Clustal parser, pass the path to the locally stored
clustal file as a parameter.

Usage:

```cpp
std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
std::size_t numberOfSequences = NUMBER_OF_SEQUENCES(myClustal);
```

---

## `ID(clustal, index)`

`ID(clustal, index)` overloaded parametric macro returns the ID of an
entry on `index`th row in clustal file guarded by `clustal`.

The parameter `clustal` is an instance of `CLUSTAL` class.

The parameter `index` is of type `std::size_t` and denotes the `index`th
entry of the clustal file guarded by `clustal`.
Indexing starts from 0.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
std::vector<std::string> ids;
for (std::size_t i = 0; i < NUMBER_OF_SEQUENCES(myClustal); ++i) {
  ids.emplace_back(ID(myClustal, i));
}
```

---

## `SEQUENCE(clustal, specifier)`

`SEQUENCE(clustal, specifier)` overloaded parametric macro returns the
sequence belonging either to `specifier`th entry of `clustal` or having
`specifier` as sequence's ID.

The parameter `clustal` is an instance of `CLUSTAL` class.

The parameter `specifier` is of type `std::size_t` and denotes the
`specifier`th entry of the clustal file guarded by `clustal`.
Indexing starts from 0.
Or is of type `std::string` and denotes the ID of the entry where returned
sequence is to be found.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
for (std::size_t i = 0; i < NUMBER_OF_SEQUENCES(myClustal); ++i) {
  assert(SEQUENCE(myClustal, i) == SEQUENCE(myClustal, ID(myClustal, i)));
}
```

---

## `COLUMN(clustal, index)`

`COLUMN(clustal, index)` parametric macro returns the `index`th column of
the multiple sequence aligment guarded by `clustal`.

The parameter `clustal` is an instance of `CLUSTAL` class.

The parameter `index` is of type `std::size_t` and denotes the `index`th
column of the MSA guarded by `clustal`.
Indexing starts from 0.

Return type is `std::vector<char>`.

Usage:

```cpp
std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
std::vector<std::vector<char>> columns;
for (std::size_t i = 0; i < NUMBER_OF_COLUMNS(myClustal); ++i) {
  columns.push_back(COLUMN(myClustal, i));
}
```

---

## `NUMBER_OF_SEQUENCES(clustal)`

`NUMBER_OF_SEQUENCES(clustal)` parametric macro returns number of entries
in the clustal file guarded by `clustal`.

The parameter `clustal` is an instance of `CLUSTAL` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
std::vector<std::string> ids;
for (std::size_t i = 0; i < NUMBER_OF_SEQUENCES(myClustal); ++i) {
  ids.emplace_back(ID(myClustal, i));
}
```

---

## `NUMBER_OF_COLUMNS(clustal)`

`NUMBER_OF_COLUMNS(clustal)` parametric macro returns number of columns
of the multiple sequence alignment guarded by `clustal`.

The parameter `clustal` is an instance of `CLUSTAL` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
std::vector<std::vector<char>> columns;
for (std::size_t i = 0; i < NUMBER_OF_COLUMNS(myClustal); ++i) {
  columns.push_back(COLUMN(myClustal, i));
}
```

---

## `SUM_OF_PAIRS(clustal, matrix, order)`

`SUM_OF_PAIRS(clustal, matrix, order)` overloaded parametric macro returns
the MSA score of the MSA guarded by `clustal`, based on a scoring matrix
`matrix` with rows and columns corresponding to residues with order given
by `order`.

The parameter `clustal` is an instance of `CLUSTAL` class.

The parameter `matrix` is of type `T<U<V>>`, where `T` and `U` are
sequence containers with defined `operator[]`, `V` is a numeric type
implicitly convertible to `double`, for example
`std::vector<std::vector<int>>`, and denotes a scoring matrix.
Its size must be 20x20 and it may be triangular.
`i`th row and `i`th column must correspond to the residue type given
in `order[i]`.
See [example](test/blosum64.hpp).

The parameter `order` is of type `T<char>`, where `T` is a sequence
container with defined `operator[]`, for example `std::vector<char>`, and
denotes the order of rows and columns of `matrix`.
Its size must be 20.
`order[i]` is the residue type, which appears in `matrix` in `i`th row
and `i`th column.
See [example](test/blosum64.hpp).

Return type is `double`.

```cpp
#include "blosum64.hpp"

std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
double sumOfPairs = SUM_OF_PAIRS(myClustal, blosum64::matrix, blosum64::order);
```

---

## `SUM_OF_PAIRS(clustal, matrix, order, index)`

`SUM_OF_PAIRS(clustal, matrix, order)` overloaded parametric macro returns
the MSA score of `index`th columns of the MSA guarded by `clustal`, based
on a scoring matrix `matrix` with rows and columns corresponding to
residues with order given by `order`.

The parameter `clustal` is an instance of `CLUSTAL` class.

The parameter `matrix` is of type `T<U<V>>`, where `T` and `U` are
sequence containers with defined `operator[]`, `V` is a numeric type
implicitly convertible to `double`, for example
`std::vector<std::vector<int>>`, and denotes a scoring matrix.
Its size must be 20x20 and it may be triangular.
`i`th row and `i`th column must correspond to the residue type given
in `order[i]`.
See [example](test/blosum64.hpp).

The parameter `order` is of type `T<char>`, where `T` is a sequence
container with defined `operator[]`, for example `std::vector<char>`, and
denotes the order of rows and columns of `matrix`.
Its size must be 20.
`order[i]` is the residue type, which appears in `matrix` in `i`th row
and `i`th column.
See [example](test/blosum64.hpp).

The parameter `index` is of type `std::size_t` and denotes the `index`th
column of `clustal`.

Return type is `double`.

```cpp
#include "blosum64.hpp"

std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
double sumOfPairsTotal = SUM_OF_PAIRS(myClustal, blosum64::matrix, blosum64::order);
double sumOfPairSum = 0;
for (std::size_t i = 0; i < NUMBER_OF_COLUMNS(myClustal); ++i) {
  sumOfPairSum += SUM_OF_PAIRS(myClustal, blosum64::matrix, blosum64::order, i);
}
assert(sumOfPairsTotal == sumOfPairSum);
```

---

## `BEST_SCORING_COLUMNS(clustal, matrix, order, n)`

`BEST_SCORING_COLUMNS(clustal, matrix, order, n)` parametric macro returns
`n` best scoring columns of the MSA guarded by `clustal`, based on a
scoring matrix `matrix` with rows and columns corresponding to residues
with order given by `order`.

The parameter `clustal` is an instance of `CLUSTAL` class.

The parameter `matrix` is of type `T<U<V>>`, where `T` and `U` are
sequence containers with defined `operator[]`, `V` is a numeric type
implicitly convertible to `double`, for example
`std::vector<std::vector<int>>`, and denotes a scoring matrix.
Its size must be 20x20 and it may be triangular.
`i`th row and `i`th column must correspond to the residue type given
in `order[i]`.
See [example](test/blosum64.hpp).

The parameter `order` is of type `T<char>`, where `T` is a sequence
container with defined `operator[]`, for example `std::vector<char>`, and
denotes the order of rows and columns of `matrix`.
Its size must be 20.
`order[i]` is the residue type, which appears in `matrix` in `i`th row
and `i`th column.
See [example](test/blosum64.hpp).

The parameter `n` is of type `std::size_t` and denotes the number of best
scoring columns to be returned.

Return type is `std::vector<std::tuple<double, std::size_t, std::shared_ptr<std::vector<char>>>>`.

Usage:

```cpp
#include "blosum64.hpp"

std::string path = "/home/hamalcij/clustal/PF00069";
CLUSTAL myClustal(path);
for (auto&& [score, index, column] : BEST_SCORING_COLUMNS(myClustal, blosum64::triangularMatrix, blosum64::order, n)) {
  std::cout << "Score: " << (int) score << "; Index: " << index;
  std::cout << "; Column: ";
  for (auto&& aminoAcid : *column) {
    std::cout << aminoAcid;
  }
  std::cout << std::endl;
}
```

---

[Previous](pdb.md)

[Up](README.md)

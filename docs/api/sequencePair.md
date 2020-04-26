[//]: # (biotool)
[//]: # (Bioinformatics Toolbox)
[//]: # ()
[//]: # (docs/api/sequencePair.md)
[//]: # (Copyright © 2020 Hamalčík Jan)
[//]: # ()
[//]: # (Documentation of SequencePair class API)
[//]: # ()

# Sequence Pair API Documentation

Following macros will conveniently let you use `biotool`'s API to work
with sequence pairs:

## `SEQUENCE_PAIR`

`SEQUENCE_PAIR` macro defines a class that will let you access methods
such as getting the alignment of two sequences using dynamic programming.
Works very well with the [Fasta API](fasta.md).

To initiate the class, pass two sequences as parameters.

Usage:

```cpp
std::string path = "/home/hamalcij/fasta/hemoglobin.fasta";
FASTA myFasta(path);
MOLECULE firstMolecule(myFasta, 0);
MOLECULE secondMolecule(myFasta, 1);
std::string sequence1 = SEQUENCE(firstMolecule);
std::string sequence2 = SEQUENCE(secondMolecule);
SEQUENCE_PAIR mySequencePair(sequence1, sequence2);
```

---

## `SEQUENCES(sequencePair)`

`SEQUENCES(sequencePair)` parametric macro returns both sequences passed
to `SEQUENCE_PAIR`'s constructor.

The parameter `sequencePair` is an instance of `SEQUENCE_PAIR` class.

Return type is `std::tuple<std::string, std::string>`.

Usage:

```cpp
std::string sequence1 = "MPRGAVILEFWTCSS";
std::string sequence2 = "MDEKECSHGGAF";
SEQUENCE_PAIR mySequencePair(sequence1, sequence2);
auto [seq1, seq2] = SEQUENCES(mySequencePair);
```

---

## `HAMMING_DISTANCE(sequencePair)`

`HAMMING_DISTANCE(sequencePair)` parametric macro returns the Hamming
distance of two sequences passed to `SEQUENCE_PAIR`'s constructor.
If the length of both sequences is not the same, the method will throw
`std::length_error`.

The parameter `sequencePair` is an instance of `SEQUENCE_PAIR` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string sequence1 = "SHORTS";
std::string sequence2 = "HORSTS";
SEQUENCE_PAIR mySequencePair(sequence1, sequence2);
std::size_t hammingDistance = HAMMING_DISTANCE(mySequencePair);
```

---

## `EDIT_DISTANCE(sequencePair)`

`EDIT_DISTANCE(sequencePair)` parametric macro returns the edit distance
of two sequences passed to `SEQUENCE_PAIR`'s constructor, as well as
all optimal alignments.
Match is scored `0`, mismatch or gapped is scored `1`.

The parameter `sequencePair` is an instance of `SEQUENCE_PAIR` class.

Return type is `std::tuple<std::size_t, std::vector<std::tuple<std::string, std::string>>>`.

Usage:

```cpp
std::string sequence1 = "SHORTS";
std::string sequence2 = "HORSTS";
SEQUENCE_PAIR mySequencePair(sequence1, sequence2);
auto [editDistance, alignments] = EDIT_DISTANCE(mySequencePair);
for (auto&& [seq1, seq2] : alignments) {
  std::cout << seq1 << std::endl;
  std::cout << seq2 << std::endl;
  std::cout << std::endl;
}
```

---

[Next](pdb.md)

[Previous](fasta.md)

[Up](README.md)

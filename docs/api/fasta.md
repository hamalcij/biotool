[//]: # (biotool)
[//]: # (Bioinformatics Toolbox)
[//]: # ()
[//]: # (docs/api/fasta.md)
[//]: # (Copyright © 2020 Hamalčík Jan)
[//]: # ()
[//]: # (Documentation of Fasta class API)
[//]: # ()

# Fasta API Documentation

Following macros will conveniently let you use `biotool`' Fasta API:

## `FASTA`

`FASTA` macro defines a class that will let you access the Fasta
parser's methods, such as getting the sequences stored in your `.fasta`
file.

To initiate the Fasta parser, pass the path to the locally stored `.fasta`
file as a parameter.

Usage:

```cpp
std::string path = "/home/hamalcij/fasta/hemoglobin.fasta";
FASTA myFasta(path);
std::size_t numberOfMolecules NUMBER_OF_MOLECULES(myFasta);
```

---

## `MOLECULE`

`MOLECULE` macro defines a class that will let you access specific
molecule's attributes from your `.fasta` file.

To instantiate the class, pass your Fasta class and index of the molecule
within the Fasta file.
Indexing starts from `0`.

Usage:

```cpp
std::string path = "/home/hamalcij/fasta/hemoglobin.fasta";
FASTA myFasta(path);
std::vector<std::string> sequences;
for (std::size_t i = 0; i < NUMBER_OF_MOLECULES(myFasta); ++i) {
  MOLECULE myMolecule(myFasta, i);
  sequences.emplace_back(SEQUENCE(myMolecule));
}
```

---

## `NUMBER_OF_MOLECULES(fasta)`

`NUMBER_OF_MOLECULES(fasta)` parametric macro returns the number of
molecules contained in the Fasta file, which is guarded by an instance of
`FASTA` class.

The parameter `fasta` is an instance of `FASTA` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/fasta/hemoglobin.fasta";
FASTA myFasta(path);
std::size_t numberOfMolecules NUMBER_OF_MOLECULES(myFasta);
```

---

## `DESCRIPTION(molecule)`

`DESCRIPTION(molecule)` parametric macro returns the description of a
molecule guarded by an instance of `MOLECULE` class.
The description corresponds to rows in the `.fasta` file starting with
`>`.

The parameter `molecule` is an instance of `MOLECULE` class.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/fasta/hemoglobin.fasta";
FASTA myFasta(path);
std::vector<std::string> descriptions;
for (std::size_t i = 0; i < NUMBER_OF_MOLECULES(myFasta); ++i) {
  MOLECULE myMolecule(myFasta, i);
  descriptions.emplace_back(DESCRIPTION(myMolecule));
}
```

---

## `SEQUENCE(molecule)`

`SEQUENCE(molecule)` overloaded parametric macro returns the sequence of
a molecule guarded by an instance of `MOLECULE` class.

The parameter `molecule` is an instance of `MOLECULE` class.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/fasta/hemoglobin.fasta";
FASTA myFasta(path);
std::vector<std::string> sequences;
for (std::size_t i = 0; i < NUMBER_OF_MOLECULES(myFasta); ++i) {
  MOLECULE myMolecule(myFasta, i);
  sequences.emplace_back(SEQUENCE(myMolecule));
}
```

---

## `SUBSEQUENCE(molecule, from, to)`

`SUBSEQUENCE(molecule, from, to)` parametric macro return the subsequence
of a molecule guarded by an instance of `MOLECULE` class starting from
residue `from` to residue `to`.

The parameter `molecule` is an instance of `MOLECULE` class.

The parameter `from` is of type `std::size_t` and denotes the `from`th
residue of `molecule`'s sequence.
`from` must be less than `to`.

The parameter `to` is of type `std::size_t` and denotes the `to`th
residue of `molecule`'s sequence.
`to` must be greater than `from`.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/fasta/hemoglobin.fasta";
FASTA myFasta(path);
std::vector<std::string> subsequences;
std::size_t from = 2;
std::size_t to = 24;
for (std::size_t i = 0; i < NUMBER_OF_MOLECULES(myFasta); ++i) {
  MOLECULE myMolecule(myFasta, i);
  subsequences.emplace_back(SUBSEQUENCE(myMolecule, from, to));
}
```

---

## `SEQUENCE_LENGTH(molecule)`

`SEQUENCE_LENGTH(molecule)` parametric macro returns the length of the
sequence of a molecule guarded by an instance of `MOLECULE` class.

The parameter `molecule` is an instance of `MOLECULE` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/fasta/hemoglobin.fasta";
FASTA myFasta(path);
std::vector<std::size_t> sequenceLengths;
for (std::size_t i = 0; i < NUMBER_OF_MOLECULES(myFasta); ++i) {
  MOLECULE myMolecule(myFasta, i);
  sequenceLengths.emplace_back(SEQUENCE_LENGTH(myMolecule));
}
```

---

[Next](sequencePair)

[Up](README.md)

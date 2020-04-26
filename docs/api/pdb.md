[//]: # (biotool)
[//]: # (Bioinformatics Toolbox)
[//]: # ()
[//]: # (docs/api/pdb.md)
[//]: # (Copyright © 2020 Hamalčík Jan)
[//]: # ()
[//]: # (Documentation of Pdb class API)
[//]: # ()

# PDB API Documentation

Following macros will conveniently let you use `biotool`' PDB API:

## `PDB`

`PDB` macro defines a class that will let you access the PDB
parser's methods, such as getting the coordinates of atom stored in your
`.pdb` file.

To initiate the PDB parser, pass the path to the locally stored `.pdb`
file as a parameter.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
std::size_t numberOfModels = NUMBER_OF_MODELS(myPdb);
```

---

## `MODEL`

`MODEL` macro defines a class that will let you access specific model's
attributes from your `.pdb` file.

To instantiate the class, use the macro `FIND_MODEL`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
```

---

## `CHAIN`

`CHAIN` macro defines a class that will let you access specific chain's
attributes from your `.pdb` file.

To instantiate the class, use the macro `FIND_CHAIN`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
```

## `RESIDUE`

`RESIDUE` macro defines a class that will let you access specific
residue's attributes from your `.pdb` file.

To instantiate the class, use the macro `FIND_RESIDUE`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
```

---

## `HET_RESIDUE`

`HET_RESIDUE` macro defines a class that will let you access specific
heterogenic residue's (ligand's) attributes from your `.pdb` file.

To instantiate the class, use the macro `FIND_HET_RESIDUE`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
```

---

## `ATOM`

`ATOM` macro defines a class that will let you access specific
atom's attributes from your `.pdb` file.

To instantiate the class, use the macro `FIND_ATOM`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOM myAtom = FIND_ATOM(myResidue, "7");
```

---

## `HET_ATOM`

`HET_ATOM` macro defines a class that will let you access specific
ligand's atom's attributes from your `.pdb` file.

To instantiate the class, use the macro `FIND_ATOM`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
```

---

## `CHAINS`

`CHAINS` macro defines a vector of `CHAIN` class instances.

To fill the vector with `CHAIN` instances, use the macro `GET_CHAINS`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAINS myChains;
GET_CHAINS(myModel, myChains);
for (auto&& chain : myChains) {
  std::cout << "Chain " << ID(chain) << " has " << NUMBER_OF_HET_RESIDUES(chain) << " ligands." << std::endl;
}
```

---

## `RESIDUES`

`RESIDUES` macro defines a vector of `RESIDUE` class instances.

To fill the vector with `RESIDUE` instances, use the macro
`GET_RESIDUES`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUES myResidues;
GET_RESIDUES(myChain, myResidues);
for (auto&& residue : myResidues) {
  std::cout << "Residue " << ID(residue) << " has " << NUMBER_OF_ATOMS(residue) << " atoms." << std::endl;
}
```

---

## `HET_RESIDUES`

`HET_RESIDUES` macro defines a vector of `HET_RESIDUE` class instances.

To fill the vector with `HET_RESIDUE` instances, use the macro
`GET_RESIDUES`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
HET_RESIDUES myLigands;
GET_HET_RESIDUES(myChain, myLigands);
for (auto&& ligand : myLigands) {
  std::cout << "Ligand " << ID(ligand) << " has " << NUMBER_OF_HET_ATOMS(ligand) << " atoms." << std::endl;
}
```

---

## `ATOMS`

`ATOMS` macro defines a vector of `ATOM` class instances.

To fill the vector with `ATOM` instances, use the macro `GET_ATOMS`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOMS myAtoms;
GET_ATOMS(myResidue, myAtoms);
for (auto&& atom : myAtoms) {
  auto [x, y, z] = COORDINATES(atom);
  std::cout << "Atom " << ID(atom) << " has following coordinates: (" << x << "," << y << ";" << z << ")" << std::endl;
}
```

---

## `HET_ATOMS`

`HET_ATOMS` macro defines a vector of `HET_ATOM` class instances.

To fill the vector with `HET_ATOM` instances, use the macro `GET_ATOMS`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "1");
HET_ATOMS myAtoms;
GET_ATOMS(myLigand, myAtoms);
for (auto&& atom : myAtoms) {
  auto [x, y, z] = COORDINATES(atom);
  std::cout << "Atom " << ID(atom) << " has following coordinates: (" << x << "," << y << ";" << z << ")" << std::endl;
}
```

---

## `FIND_MODEL(pdb, id)`

`FIND_MODEL(pdb, id)` parametric macro returns an instance of `MODEL`
class from the `.pdb` file guarded by an instance of `PDB` class.

The parameter `pdb` is an instance of `PDB` class.

The parameter `id` is of type `std::string` and denotes the model's
unique ID.

Return type is `MODEL`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
std::size_t numberOfChains = NUMBER_OF_CHAINS(myModel);
```

---

## `FIND_CHAIN(model, id)`

`FIND_CHAIN(model, id)` parametric macro returns an instance of `CHAIN`
class from the `.pdb` file guarded by an instance of `PDB` class.

The parameter `model` is an instance of `MODEL` class.

The parameter `id` is of type `std::string` and denotes the chain's
unique ID.

Return type is `CHAIN`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
std::size_t numberOfLigands = NUMBER_OF_HET_RESIDUES(myChain);
```

---

## `FIND_RESIDUE(chain, id)`

`FIND_RESIDUE(chain, id)` parametric macro returns an instance of `RESIDUE`
class from the `.pdb` file guarded by an instance of `PDB` class.

The parameter `chain` is an instance of `CHAIN` class.

The parameter `id` is of type `std::string` and denotes the residue's
unique ID.

Return type is `RESIDUE`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
std::size_t numberOfAtoms = NUMBER_OF_ATOMS(myResidue);
```

---

## `FIND_HET_RESIDUE(chain, id)`

`FIND_HET_RESIDUE(chain, id)` parametric macro returns an instance of
`HET_RESIDUE`
class from the `.pdb` file guarded by an instance of `PDB` class.

The parameter `chain` is an instance of `CHAIN` class.

The parameter `id` is of type `std::string` and denotes the ligand's
unique ID.

Return type is `HET_RESIDUE`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
std::size_t numberOfAtoms = NUMBER_OF_ATOMS(myLigand);
```

---

## `FIND_ATOM(residue, id)`

`FIND_ATOM(residue, id)` overloaded parametric macro returns an instance of `ATOM` or `HET_ATOM`
class from the `.pdb` file guarded by an instance of `PDB` class.

The parameter `residue` is an instance of `RESIDUE` or `HET_RESIDUE`
class.

The parameter `id` is of type `std::string` and denotes the atom's
unique ID.

Return type is `ATOM` or `HET_ATOM`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
ATOM myAtom = FIND_ATOM(myResidue, "7");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
auto [x1, y1, z1] = COORDINATES(myAtom);
auto [x2, y2, z2] = COORDINATES(myLigandAtom);
```

---

## `GET_MODELS(pdb)`

`GET_MODELS(pdb)` parametric macro returns a vector of `MODEL` instances
from the `.pdb` file guarded by an instance of `PDB` class.

The parameter `pdb` is an instance of `PDB` class.

Return type is `std::vector<MODEL>`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
for (auto&& model : GET_MODELS(myPdb)) {
  std::cout << "Model " << ID(model) << " has " << NUMBER_OF_CHAINS(model) << " chains." << std::endl;
}
```

---

## `GET_CHAINS(model, chains)`

`GET_CHAINS(model, chains)` parametric macro emplaces all chains of
`model` into `chains`.

The parameter `model` is an instance of `MODEL` class.

The parameter `chains` is an instance of `CHAINS` vector.
Any data which it contains will be deleted before generating new data.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
for (auto&& model : GET_MODELS(myPdb)) {
  CHAINS myChains;
  GET_CHAINS(model, myChains);
  for (auto&& chain : myChains) {
    std::cout << "Chain " << ID(chain) << " has " << NUMBER_OF_RESIDUES(chain) << " residues." << std::endl;
  }
}
```

---

## `GET_RESIDUES(chain, residues)`

`GET_RESIDUES(chain, residues)` parametric macro emplaces all residues of
`chain` into `residues`.

The parameter `chain` is an instance of `CHAIN` class.

The parameter `residues` is an instance of `RESIDUES` vector.
Any data which it contains will be deleted before generating new data.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
for (auto&& model : GET_MODELS(myPdb)) {
  CHAINS myChains;
  GET_CHAINS(model, myChains);
  for (auto&& chain : myChains) {
    RESIDUES myResidues;
    GET_RESIDUES(chain, myResidues);
    for (auto&& residue : myResidues) {
      std::cout << "Residue " << ID(residue) << " has " << NUMBER_OF_ATOMS(residue) << " atoms." << std::endl;
    }
  }
}
```

---

## `GET_HET_RESIDUES(chain, hetResidues)`

`GET_HET_RESIDUES(chain, hetResidues)` parametric macro emplaces all
heterogenic residues (ligands) of `chain` into `hetResidues`.

The parameter `chain` is an instance of `CHAIN` class.

The parameter `hetResidues` is an instance of `HET_RESIDUES` vector.
Any data which it contains will be deleted before generating new data.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
for (auto&& model : GET_MODELS(myPdb)) {
  CHAINS myChains;
  GET_CHAINS(model, myChains);
  for (auto&& chain : myChains) {
    HET_RESIDUES myLigands;
    GET_HET_RESIDUES(chain, myLigands);
    for (auto&& ligand : myLigands) {
      std::cout << "Ligand " << ID(ligand) << " has " << NUMBER_OF_HET_ATOMS(ligand) << " atoms." << std::endl;
    }
  }
}
```

---

## `GET_ATOMS(residue, atoms)`

`GET_ATOMS(residue, atoms)` overloaded parametric macro emplaces all
atoms of `residue` into `atoms`.

The parameter `residue` is an instance of `RESIDUE` or `HET_RESIDUE`
class.

The parameter `atoms` is an instance of `ATOMS` or `HET_ATOMS` vector.
Any data which it contains will be deleted before generating new data.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
for (auto&& model : GET_MODELS(myPdb)) {
  CHAINS myChains;
  GET_CHAINS(model, myChains);
  for (auto&& chain : myChains) {
    RESIDUES myResidues;
    GET_RESIDUES(chain, myResidues);
    for (auto&& residue : myResidues) {
      ATOMS myAtoms;
      GET_ATOMS(residue, myAtoms);
      for (auto&& atom : myAtoms) {
        auto [x, y, z] = COORDINATES(atom);
        std::cout << "Atom " << ID(atom) << " has following coordinates: (" << x << "," << y << ";" << z << ")" << std::endl;
      }
    }
    HET_RESIDUES myLigands;
    GET_HET_RESIDUES(chain, myLigands);
    for (auto&& ligand : myLigands) {
      HET_ATOMS myLigandAtoms;
      GET_ATOMS(ligand, myLigandAtoms);
      for (auto&& atom : myAtoms) {
        auto [x, y, z] = COORDINATES(atom);
        std::cout << "Ligand atom " << ID(atom) << " has following coordinates: (" << x << "," << y << ";" << z << ")" << std::endl;
      }
    }
  }
}
```

---

## `NUMBER_OF_MODELS(pdb)`

`NUMBER_OF_MODELS(pdb)` parametric macro returns the number of models
contained in a `.pdb` file guarded by class `PDB`.

The parameter `pdb` is an instance of `PDB` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
std::size_t numberOfModels = NUMBER_OF_MODELS(myPdb);
```

---

## `NUMBER_OF_CHAINS(model)`

`NUMBER_OF_CHAINS(model)` parametric macro returns the number of chains
contained in `model`.

The parameter `model` is an instance of `MODEL` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
std::size_t numberOfChains = NUMBER_OF_CHAINS(myModel);
```

---

## `NUMBER_OF_RESIDUES(unit)`

`NUMBER_OF_RESIDUES(unit)` overloaded parametric macro returns the number
of residues contained in `unit`.

The parameter `unit` is an instance of `MODEL` or `CHAIN` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
std::size_t numberOfResiduesInModel = NUMBER_OF_RESIDUES(myModel);
std::size_t numberOfResiduesInChain = NUMBER_OF_RESIDUES(myChain);
```

---

## `NUMBER_OF_HET_RESIDUES(unit)`

`NUMBER_OF_HET_RESIDUES(unit)` overloaded parametric macro returns the
number of heterogenic residues (ligands) contained in `unit`.

The parameter `unit` is an instance of `MODEL` or `CHAIN` class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
std::size_t numberOfLigandsInModel = NUMBER_OF_HET_RESIDUES(myModel);
std::size_t numberOfLigandsInChain = NUMBER_OF_HET_RESIDUES(myChain);
```

---

## `NUMBER_OF_ATOMS(unit)`

`NUMBER_OF_ATOMS(unit)` overloaded parametric macro returns the
number of atoms contained in `unit`.

The parameter `unit` is an instance of `MODEL` or `CHAIN` or `RESIDUE`
class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
std::size_t numberOfAtomsInModel = NUMBER_OF_ATOMS(myModel);
std::size_t numberOfAtomsInChain = NUMBER_OF_ATOMS(myChain);
std::size_t numberOfAtomsInResidue = NUMBER_OF_ATOMS(myResidue);
```

---

## `NUMBER_OF_HET_ATOMS(unit)`

`NUMBER_OF_HET_ATOMS(unit)` overloaded parametric macro returns the
number of heterogenic residue's (ligand's) atoms contained in `unit`.

The parameter `unit` is an instance of `MODEL` or `CHAIN` or `HET_RESIDUE`
class.

Return type is `std::size_t`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
std::size_t numberOfLigandAtomsInModel = NUMBER_OF_HET_ATOMS(myModel);
std::size_t numberOfLigandAtomsInChain = NUMBER_OF_HET_ATOMS(myChain);
std::size_t numberOfLigandAtomsInResidue = NUMBER_OF_HET_ATOMS(myResidue);
```

---

## `ID(unit)`

`ID(unit)` overloaded parametric macro returns the ID of `unit`.

The parameter `unit` is an instance of `MODEL` or `CHAIN` or `RESIDUE`
or `HET_RESIDUE` or `ATOM` or `HET_ATOM` class.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOM myAtom = FIND_ATOM(myResidue, "7");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
std::string modelID = ID(myModel);
std::string chainID = ID(myChain);
std::string residueID = ID(myResidue);
std::string atomID = ID(myAtom);
std::string ligandID = ID(myLigand);
std::string ligandAtomID = ID(myLigandAtom);
```

---

## `RESIDUE_NAME(residue)`

`RESIDUE_NAME(residue)` overloaded parametric macro returns the name of `residue`.

The parameter `residue` is an instance of `RESIDUE` or `HET_RESIDUE`
class.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
std::string residueName = RESIDUE_NAME(myResidue);
std::string ligandName = RESIDUE_NAME(myLigand);
```

---

## `ATOM_NAME(atom)`

`ATOM_NAME(atom)` overloaded parametric macro returns the name of `atom`.

The parameter `atom` is an instance of `ATOM` or `HET_ATOM` class.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOM myAtom = FIND_ATOM(myResidue, "7");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
std::string atomName = ATOM_NAME(myAtom);
std::string ligandAtomName = ATOM_NAME(myLigandAtom);
```

---

## `COORDINATES(atom)`

`COORDINATES(atom)` overloaded parametric macro returns the coordinates of
`atom`.

The parameter `atom` is an instance of `ATOM` or `HET_ATOM` class.

Return type is `std::tuple<float, float, float>`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOM myAtom = FIND_ATOM(myResidue, "7");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
auto [x1, y1, z1] atomCoords = COORDINATES(myAtom);
auto [x2, y2, z2] ligandAtomCoords = COORDINATES(myLigandAtom);
```

---

## `OCCUPANCY(atom)`

`OCCUPANCY(atom)` overloaded parametric macro returns the occupancy of
`atom`.

The parameter `atom` is an instance of `ATOM` or `HET_ATOM` class.

Return type is `float`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOM myAtom = FIND_ATOM(myResidue, "7");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
float atomOccupancy = OCCUPANCY(myAtom);
float ligandAtomOccupancy = OCCUPANCY(myLigandAtom);
```

---

## `TEMPERATURE_FACTOR(atom)`

`TEMPERATURE_FACTOR(atom)` overloaded parametric macro returns the
temperatur factor of `atom`.

The parameter `atom` is an instance of `ATOM` or `HET_ATOM` class.

Return type is `float`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOM myAtom = FIND_ATOM(myResidue, "7");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
float atomTempFactor = TEMPERATURE_FACTOR(myAtom);
float ligandAtomTempFactor = TEMPERATURE_FACTOR(myLigandAtom);
```

---

## `ELEMENT(atom)`

`ELEMENT(atom)` overloaded parametric macro returns the element of `atom`.

The parameter `atom` is an instance of `ATOM` or `HET_ATOM` class.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOM myAtom = FIND_ATOM(myResidue, "7");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
std::string atomElement = ELEMENT(myAtom);
std::string ligandAtomElement = ELEMENT(myLigandAtom);
```

---

## `CHARGE(atom)`

`CHARGE(atom)` overloaded parametric macro returns the charge of `atom`.

The parameter `atom` is an instance of `ATOM` or `HET_ATOM` class.

Return type is `std::string`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
RESIDUE myResidue = FIND_RESIDUE(myChain, "1");
ATOM myAtom = FIND_ATOM(myResidue, "7");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
HET_ATOM myLigandAtom = FIND_ATOM(myLigand, "4438");
std::string atomCharge = CHARGE(myAtom);
std::string ligandAtomCharge = CHARGE(myLigandAtom);
```

---

## `GET_ATOMS_CLOSE_TO_LIGAND(model, ligand, atoms, maxDistance)`

`GET_ATOMS_CLOSE_TO_LIGAND(model, ligand, atoms, maxDistance)` parametric
macro assigns `ATOM` instances that are within `maxDistance` from
`HET_ATOM` instances belonging to `ligand` in `model`.

The parameter `model` is an instance of `MODEL` class.

The parameter `ligand` is an instance of `HET_RESIDUE` class.

The parameter `atoms` is an instance of `ATOMS` vector.
Any data which it contains will be deleted before generating new data.

The parameter `maxDistance` is of type `float` and denotes the maximum
distance in Ånströms from `ATOM` instances belonging to `ligand` to any
`ATOM` instance which should be inserted to `atoms`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
ATOMS nearAtoms;
float maxDistance = 3.2;
GET_ATOMS_CLOSE_TO_LIGAND(myModel, myLigand, nearAtoms, maxDistance);
for (auto&& atom : nearAtoms) {
  auto [x, y, z] = COORDINATES(atom);
  std::cout << "Atom " << ID(atom) << " has following coordinates: (" << x << "," << y << ";" << z << ")" << std::endl;
}
```

---

## `GET_RESIDUES_CLOSE_TO_LIGAND(model, ligand, residues, maxDistance)`

`GET_RESIDUES_CLOSE_TO_LIGAND(model, ligand, residues, maxDistance)`
parametric macro assigns `RESIDUE` instances that are within `maxDistance`
from `HET_ATOM` instances belonging to `ligand` in `model`.

The parameter `model` is an instance of `MODEL` class.

The parameter `ligand` is an instance of `HET_RESIDUE` class.

The parameter `residues` is an instance of `RESIDUES` vector.
Any data which it contains will be deleted before generating new data.

The parameter `maxDistance` is of type `float` and denotes the maximum
distance in Ånströms from `ATOM` instances belonging to `ligand` to any
`ATOM` instance belonging to a `RESIDUE` instance which should be
inserted to `atoms`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
CHAIN myChain = FIND_CHAIN(myModel, "A");
HET_RESIDUE myLigand = FIND_HET_RESIDUE(myChain, "156");
RESIDUES nearResidues;
float maxDistance = 3.2;
GET_RESIDUES_CLOSE_TO_LIGAND(myModel, myLigand, nearResidues, maxDistance);
for (auto&& residue : nearResidues) {
  std::cout << "Residue " << ID(residue) << " has " << NUMBER_OF_ATOMS(residue) << " atoms." << std::endl;
}
```

---

## `NUMBER_OF_SURFACE_AND_BURIED_RESIDUES(model)`

`NUMBER_OF_SURFACE_AND_BURIED_RESIDUES(model)` parametric macro returns
the number of residues of `model` that have at least one atom accessible
to the solvent (a water molecule with diameter 2.75Å) as well as the
number of residues of `model` that are buried within the core of the
protein.
Notice that this is only an approximate method!

The parameter `model` is an instance of `MODEL` class.

Return type is `std::tuple<std::size_t, std::size_t>`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
auto [surface, buried] = NUMBER_OF_SURFACE_AND_BURIED_RESIDUES(myModel);
```

---

## `GET_SURFACE_AND_BURIED_STATS(model)`

`GET_SURFACE_AND_BURIED_STATS(model)` parametric macro returns the
statistic of how many of each 20 different residue types have at least
one atom accessible to the solvent (a water molecule with diameter 2.75Å)
is `model` as well as how many of each 20 different residue types are
buried within the core of `model`.
Notice that this is only an approximate method!

The parameter `model` is an instance of `MODEL` class.

Return type is `std::tuple<std::unordered_map<std::string, std::size_t>, std::unordered_map<std::string, std::size_t>>`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
auto [surfaceStats, buriedStats] = GET_PORTION_OF_POLAR_SURFACE_AND_BURIED(myModel);
for (auto&& [name, count] : surfaceStats) {
  std::cout << name << " " << count << std::endl;
}
for (auto&& [name, count] : buriedStats) {
  std::cout << name << " " << count << std::endl;
}
```

---

## `GET_PORTION_OF_POLAR_SURFACE_AND_BURIED(model)`

`GET_PORTION_OF_POLAR_SURFACE_AND_BURIED(model)` parametric macro returns
the portion of polar residues among the residues that have at least one
atom accessible to the solvent (a water molecule with diameter 2.75Å) as
well as the portion of polar residues among the residues that are buried
within the core of `model`.
Notice that this is only an approximate method!

The parameter `model` is an instance of `MODEL` class.

Return type is `std::tuple<float, float>`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
auto [polarSurfacePortion, polarBuriedPortion] = GET_PORTION_OF_POLAR_SURFACE_AND_BURIED(myModel);
```

---

## `MODEL_WIDTH(model)`

`MODEL_WIDTH(model)` parametric macro returns the distance between two
farthest atoms of `model`.

The parameter `model` is an instance of `MODEL` class.

Return type is `float`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
float width = MODEL_WIDTH(myModel);
```

---

## `MODEL_DIAMETER(model)`

`MODEL_DIAMETER(model)` parametric macro return the diameter of the
smallest circumsphere of `model`'s atoms.
Notice that you may need to run this method more than once and fetch the
lowest result as there exists more than one circumsphere which may seem
minimal.

The parameter `model` is an instance of `MODEL` class.

Return type is `float`.

Usage:

```cpp
std::string path = "/home/hamalcij/pdb/1B0B.pdb";
PDB myPdb(path);
MODEL myModel = FIND_MODEL(myPdb, "1");
float diameter = MODEL_DIAMETER(myModel);
```

---

[Next](clustal.md)

[Previous](sequencePair.md)

[Up](README.md)

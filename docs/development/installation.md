[//]: # (biotool)
[//]: # (Bioinformatics Toolbox)
[//]: # ()
[//]: # (docs/development/installation.md)
[//]: # (Copyright © 2020 Hamalčík Jan)
[//]: # ()
[//]: # (Informs about how to use biotool in your project)
[//]: # ()

# Installation

If the requirements specified [here](requirements.md) are being met for
your project, you may continue with the installation of `biotool`.
First of all, download the
[archive](https://github.com/hamalcij/biotool/archive/master.zip) with
`biotool`'s source code from GitHub, navigate to your Downloads folder,
extract the `.zip` file and move the extracted directory to the local
`git` repository with your C++ bioinformatics project by typing into the
command line:

```shell
unzip biotool-master.zip
mv biotool-master /path/to/your/project/biotool
cd /path/to/your/project/biotool
```

From now on, I expect that you are in your `/path/to/your/project/biotool`
folder.

You may run and inspect the tests in `test` to see how `biotool`
C++ API may be used.
Compile them by typing `make`, which will build the tests with
instructions from `Makefile`.
If you no longer need the tests, type `rm -r test Makefile`.
You can also remove the documentation by typing `rm -r docs` to save some
space, but I do not recommend it.

## How to include `biotool`

To include `biotool` C++ API in your project, you only need to include
the `biotool.hpp` file, which in turn includes all the other necessary
source files.

When building your project, it is necessary to provide paths to all
`biotool`'s `.cpp` files in your `Makefile`.
Add the following paths to your source file list (e.g. `SRC` variable):

```Makefile
biotool/library/quickhull/QuickHull.cpp $(wildcard library/*.cpp)
```

---

[Previous](requirements.md)

[Up](README.md)

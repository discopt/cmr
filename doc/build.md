# Build instructions {#build}

The library can be built using the [CMake build system](https://cmake.org/).
A simple command sequence on Linux looks like this:

    mkdir build
    cd build
    cmake ..
    make

Below you can find a list of all parameters for building CMR.
As an optional step, one can execute unittests:

    make test

The executables now reside in the current directory. One can also install them via

    make install

**Build Parameters**

  - `-DCMAKE_BUILD_TYPE=Release` Compiles the code with optimization turned on. Make sue to use this if you need fast code for large matrices.
  - `-DGENERATORS=on`            Builds [generator tools](\ref generators) for certain matrices.
  - `-DGMP=off`                  Disables large numbers; see \ref equimodular.


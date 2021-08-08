# Combinatorial Matrix Recognition Library # {#index}

This C/C++ library contains efficient implementations of various algorithm for recognizing combinatorial matrices.
The following matrix classes can be recognized:

  - \ref TU
  - \ref NETWORK

Moreover, base representation matrices for the following matroid classes can be recognized:

  - \ref REGULAR
  - \ref GRAPHIC
  - \ref SP

The following matrix/matroid classes are planned:

  - \ref BALANCED
  - \ref PERFECT

# Usage #

The functionality can be used via command-line tools or via its C/C++ interface.
There is one executable per matrix/matroid class.
Its documentation and that of the interface can be found on the respective pages.
The executables accept several \ref FILEFORMATS.

# License #

The software is released under the [MIT License](https://en.wikipedia.org/wiki/MIT_License).
For recognition of \ref TU, please cite the paper that describes the algorithm and its implementation:

    @Article{WalterT13,
      author    = {Walter, Matthias and Truemper, Klaus},
      title     = {{Implementation of a unimodularity test}},
      journal   = {Mathematical Programming Computation},
      year      = {2013},
      volume    = {5},
      number    = {1},
      pages     = {57--73},
      issn      = {1867-2949},
      doi       = {10.1007/s12532-012-0048-x},
      publisher = {Springer-Verlag},
    }

# Authors #

- [Matthias Walter](https://people.utwente.nl/m.walter) (maintainer)
- [Klaus Truemper](https://personal.utdallas.edu/~klaus/)


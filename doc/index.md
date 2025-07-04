# Combinatorial Matrix Recognition Library # {#index}

This C/C++ library contains efficient implementations of various algorithms for recognizing combinatorial matrices.
It is the successor of the [TUtest library](\ref TUtest).
The following matrix classes can be recognized:

  - \ref tu
  - \ref network
  - \ref ctu
  - \ref equimodular
  - \ref balanced (currently only via enumerative methods)

Moreover, [representation matrices](\ref matroids) for the following matroid classes can be recognized:

  - [Binary](\ref binary_regular) and [ternary](\ref tu) regular matroids
  - [Graphic, cographic and planar](\ref graphic) matroids
  - [Series-parallel](\ref series-parallel) matroids

The following additional functionality is also part of the library:

  - \ref utilities
  - \ref k_ary
  - \ref camion
  - \ref named

The library also contains various [instance generators](\ref generators).

The following matrix/matroid classes are **planned**:

  - \ref totally-balanced
  - \ref perfect
  - \ref consecutive-ones
  - \ref max-flow-min-cut

# Installation and Usage #

The library can be found on [github](https://github.com/discopt/cmr/), or directly as a [ZIP Archive](https://github.com/discopt/cmr/archive/refs/heads/master.zip).
The functionality can be used via command-line tools or via its C interface (see [build instructions](\ref build) for more information).
There is one executable per matrix/matroid class.
Its documentation and that of the interface can be found on the respective pages.
The executables accept several \ref file-formats.

Changes and bug fixes are listed [here](\ref changes).

# License #

The software is released under the [MIT License](https://en.wikipedia.org/wiki/MIT_License).
Please cite the paper(s) corresponding to the respective tools.

# Authors # {#authors}

- [Rolf van der Hulst](https://people.utwente.nl/r.p.vanderhulst) (graphicness, network matrices)
- Henk Kraaij (balancedness)
- [Klaus Truemper](https://personal.utdallas.edu/~klaus/) (total unimodularity, regularity, complement total unimodularity)
- [Matthias Walter](https://people.utwente.nl/m.walter) (maintainer)

# Funding support 

- The software library acknowledges support from the Dutch Research Council (NWO) on grant number OCENW.M20.151 (11/2022-10/2026).

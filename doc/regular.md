# Regular Matroids # {#regular}

A matrix \f$ M \in \{0,1\}^{m \times n} \f$ is **regular** if the binary matroid [represented](\ref matroids) by \f$ M \f$ is representable over *every* field.
This is equivalent to requiring that the matroid is equal to the ternary matroid [represented](\ref matroids) by the [Camion-signed](\ref camion) version \f$ M' \f$ of \f$ M \f$, and thus equivalent to [total unimodularity](\ref tu) of \f$ M' \f$.


## Usage ##

The executable `cmr-regular` determines whether the support matrix \f$ M \f$ of a given matrix is regular.

    ./cmr-regular [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output matrices; default: `dense`.
  - `-d`        Output the decomposition tree if \f$ M \f$ is regular.
  - `n`         Output the elements of a minimal non-regular submatrix.
  - `N`         Output a minimal non-regular submatrix.
  - `s`         Print statistics about the computation to stderr.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.

## Algorithm ##

The implemented recognition algorithm is based on [Implementation of a unimodularity test](https://doi.org/10.1007/s12532-012-0048-x) by Matthias Walter and Klaus Truemper (Mathematical Programming Computation, 2013).
It is based on Seymour's [decomposition theorem for regular matroids](https://doi.org/10.1016/0095-8956(80)90075-1).
The algorithm runs in \f$ \mathcal{O}( (m+n)^5 ) \f$ time and is a simplified version of [Truemper's cubic algorithm](https://doi.org/10.1016/0095-8956(90)90030-4).
Please cite the following paper in case the implementation contributed to your research:

    @Article{WalterT13,
      author    = {Walter, Matthias and Truemper, Klaus},
      title     = {Implementation of a unimodularity test},
      journal   = {Mathematical Programming Computation},
      year      = {2013},
      volume    = {5},
      number    = {1},
      pages     = {57--73},
      issn      = {1867-2949},
      doi       = {10.1007/s12532-012-0048-x},
      publisher = {Springer-Verlag},
    }

## C Interface ##

The functionality is defined in \ref regular.h.
The main functions are:

  - CMRtestBinaryRegular() tests a binary matrix for regularity.

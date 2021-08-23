# Regular Matroids # {#regular}

A matrix \f$ M \in \{0,1\}^{m \times n} \f$ is **regular** if the binary matroid represented by \f$ M \f$ is representable over *every* field.
This is equivalent to requiring that the matroid is equal to the ternary matroid represented by the [Camion-signed](\ref camion) version \f$ M' \f$ of \f$ M \f$.


## Usage ##

The executable `cmr-regular` determines whether the support matrix \f$ M \f$ of a given matrix is regular.

    ./cmr-regular [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output matrices; default: `dense`.
  - `-d`        Output the decomposition tree if \f$ M \f$ is regular.
  - `-s`        Output the elements of a minimal non-regular submatrix.
  - `-S`        Output a minimal non-regular submatrix.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.

## Algorithm ##

The implemented recognition algorithm is based on [Implementation of a unimodularity test](https://doi.org/10.1007/s12532-012-0048-x) by Matthias Walter and Klaus Truemper (Mathematical Programming Computation, 2013).
Please cite the paper in case the implementation contributed to your research:

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

  - CMRtestRegular() tests a matrix for having a regular support matrix.

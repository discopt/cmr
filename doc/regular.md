# Regular Matroids # {#regular}

A matrix \f$ M \in \{0,1\}^{m \times n} \f$ is **regular** if the binary matroid [represented](\ref matroids) by \f$ M \f$ is representable over *every* field.
This is equivalent to requiring that the matroid is equal to the ternary matroid [represented](\ref matroids) by the [Camion-signed](\ref camion) version \f$ M' \f$ of \f$ M \f$, and thus equivalent to [total unimodularity](\ref tu) of \f$ M' \f$.


## Recognizing Regular Matroids ##

The command

    cmr-regular IN-MAT [OPTION...]

determines whether the matrix given in file `IN-MAT` is regular.

**Options:**
  - `-i FORMAT`    Format of file `IN-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-D OUT-DEC`   Write a decomposition tree of the regular matroid to file `OUT-DEC`; default: skip computation.
  - `-N NON-MINOR` Write a minimal non-regular submatrix to file `NON-SUB`; default: skip computation.
  - `-s`           Print statistics about the computation to stderr.

**Advanced options:**
  - `--no-direct-graphic`  Check only 3-connected matrices for regularity.
  - `--no-series-parallel` Do not allow series-parallel operations in decomposition tree.

If `IN-MAT` is `-` then the matrix is read from stdin.
If `OUT-DEC` or `NON-SUB` is `-` then the decomposition tree (resp. the submatrix) is written to stdout.

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

The corresponding function in the library is

  - CMRtestBinaryRegular() tests a binary matrix for regularity.

and is defined in \ref regular.h.

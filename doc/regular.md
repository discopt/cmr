# Binary Regular Matrices # {#binary_regular}

A matrix \f$ M \in \{0,1\}^{m \times n} \f$ is **regular** if the binary matroid [represented](\ref matroids) by \f$ M \f$ is representable over *every* field.
This is equivalent to requiring that the matroid is equal to the ternary matroid [represented](\ref matroids) by the [Camion-signed](\ref camion) version \f$ M' \f$ of \f$ M \f$, and thus equivalent to [total unimodularity](\ref tu) of \f$ M' \f$.


## Recognizing Binary Regular Matrices ##

The command

    cmr-regular IN-MAT [OPTION...]

determines whether the [matrix](\ref file-formats-matrix) given in file `IN-MAT` is regular.

**Options:**
  - `-i FORMAT`    Format of file `IN-MAT`; default: [dense](\ref dense-matrix).
  - `-D OUT-DEC`   Write a decomposition tree of the regular matroid to file `OUT-DEC`; default: skip computation.
  - `-N NON-MINOR` Write a minimal non-regular submatrix to file `NON-SUB`; default: skip computation.

**Advanced options:**
  - `--stats`              Print statistics about the computation to stderr.
  - `--time-limit LIMIT`   Allow at most `LIMIT` seconds for the computation.
  - `--decompose STRATEGY` Strategy for decomposing among {`DP`, `YP`, `P3`, `D3`, `Y3`}; default: `D3`.
  - `--no-direct-graphic`  Check only 3-connected matrices for regularity.
  - `--no-series-parallel` Do not allow series-parallel operations in decomposition tree.
  - `--no-simple-3-sepa`   Do not allow testing for simple 3-separations.


**Decomposition strategies:** 1st letter for distributed, 2nd for concentrated rank(s).
  - `D` Delta-sum (distributed ranks)
  - `Y` Y-sum (distributed ranks)
  - `3` 3-sum (concentrated rank)
  - `P` pivot (changes rank type)
Note that D3 and Y3 do not produce pivots.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

If `IN-MAT` is `-` then the [matrix](\ref file-formats-matrix) is read from stdin.

If `OUT-DEC` or `NON-SUB` is `-` then the decomposition tree (resp. the [submatrix](\ref file-formats-submatrix)) is written to stdout.

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

  - CMRregularTest() tests a binary matrix for regularity.

and is defined in \ref regular.h.

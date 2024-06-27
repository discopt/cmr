# Balanced / Balanceable Matrices # {#balanced}

A ternary matrix \f$ M \in \{-1,0,+1\}^{m \times n} \f$ is called **balanced** if it does not contain a square submatrix with two nonzero entries per row and per column in which the sum of all entries is 2 modulo 4.
A binary matrix \f$ M \in \{0,1\}^{m \times n} \f$ is called **balanceable** if and its nonzero entries can be signed so that the resulting matrix is balanced.

## Recognizing Balanced Matrices ##

The command

    cmr-balanced IN-MAT [OPTION...]

determines whether the [matrix](\ref file-formats-matrix) given in file `IN-MAT` is balanced.

**Options:**
  - `-i FORMAT`    Format of file `IN-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-N NON-SUB`   Write a minimal non-balanced submatrix to file `NON-SUB`; default: skip computation.
  - `-s`           Print statistics about the computation to stderr.

**Advanced options:**
  - `--time-limit LIMIT` Allow at most LIMIT seconds for the computation.
  - `--algorithm ALGO`   Algorithm to use, among `enumerate` and `graph`; default: choose best.

If `IN-MAT` is `-` then the [matrix](\ref file-formats-matrix) is read from stdin.
If `NON-SUB` is `-` then the [submatrix](\ref file-formats-submatrix) is written to stdout.

## Algorithm ##

Two types of algorithms exist for both, the binary and the ternary case.
The first algorithm types enumerate subsets of rows and then finds subsets of columns such that the resulting submatrix defines an odd cycle.
Its running time is exponential in the size of the matrix.
The second algorithm types are the [polynomial-time algorithms](https://doi.org/10.1016/j.jctb.2005.02.006) by Giacomo Zambelli (Journal of Combinatorial Theory, Series B, 2005).
They run in \f$ \mathcal{O}( (m+n)^9 ) \f$ time for binary matrices and in \f$ \mathcal{O}( (m+n)^{11} ) \f$ time for ternary matrices.

\note Only the first algorithm is implemented so far.

## C Interface ##

The corresponding function in the library is

  - CMRbalancedTest() tests a matrix for balancedness.

and is defined in \ref balanced.h.
  

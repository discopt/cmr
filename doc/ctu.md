# Complement Totally Unimodular Matrices # {#ctu}

Consider a binary matrix \f$ M \in \{0,1\}^{m \times n} \f$.
A **row complement** for row \f$ i \f$ is obtained by complementing all entries \f$ M_{r,c} \f$ with \f$ r \neq i \f$
and \f$ M_{i,c} = 1 \f$.
A **column complement** is defined as a row complement on the transpose matrix.
A binary matrix \f$ M \f$ is **complement totally unimodular** if all matrices obtainable by successive row- and column complements are [totally unimodular](\ref tu).
It turns out that all such matrices can be obtained already by at most one row complement and at most one column complement.
Hence, complement total unimodularity can be checked by checking \f$ (m+1) \cdot (n+1) \f$ matrices for [total unimodularity](\ref tu).

## Usage ##

The executable `cmr-ctu` determines whether a given matrix \f$ M \f$ is complement totally unimodular or applies
row- or column-complement operations (or both at the same time) to \f$ M \f$.

    ./cmr-ctu [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output matrices; default: `dense`.
  - `-r ROW`    Perform a row complement operation on \f$ M \f$ and do not test for complement total unimodularity.
  - `-c COLUMN` Perform a row complement operation on \f$ M \f$ and do not test for complement total unimodularity.
  - `-b`        Output the complement operations that leads to a non-totally-unimodular matrix (if \f$ M \f$ is not complement totally unimodular).
  - `-B`        Output the complemented matrix that is non-totally-unimodular (if \f$ M \f$ is not complement totally unimodular).

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.
If `-r` or `-c` (or both) are specified, then \f$ M \f$ is not tested for complement total unimodularity.

## C Interface ##

The functionality is defined in \ref ctu.h.
The main functions are:

  - CMRcomplementRowColumn() carries out a row- and column-complement operations for a matrix.
  - CMRtestComplementTotalUnimodularity() tests a matrix for being complement totally unimodular.

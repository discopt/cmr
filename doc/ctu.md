# Complement Totally Unimodular Matrices # {#ctu}

Consider a binary matrix \f$ M \in \{0,1\}^{m \times n} \f$.
A **row complement** for row \f$ i \f$ is obtained by complementing all entries \f$ M_{r,c} \f$ with \f$ r \neq i \f$
and \f$ M_{i,c} = 1 \f$.
A **column complement** is defined as a row complement on the transpose matrix.
A binary matrix \f$ M \f$ is **complement totally unimodular** if all matrices obtainable by successive row- and column complements are [totally unimodular](\ref tu).
It turns out that all such matrices can be obtained already by at most one row complement and at most one column complement.
Hence, complement total unimodularity can be checked by checking \f$ (m+1) \cdot (n+1) \f$ matrices for [total unimodularity](\ref tu).

## Recognizing Complement Totally Unimodular Matrices ##

The command

    cmr-ctu IN-MAT [OPTION]...

determines whether the [matrix](\ref file-formats-matrix) given in file `IN-MAT` is complement totally unimodular.

**Options**:
  - `-i FORMAT`   Format of file `IN-MAT`; default: [dense](\ref dense-matrix).
  - `-o FORMAT`   Format of file `OUT-MAT`; default: same as for `IN-MAT`.
  - `-n OUT-OPS`  Write complement operations that leads to a non-totally-unimodular matrix to file `OUT-OPS`; default: skip computation.
  - `-N OUT-MAT`  Write a complemented matrix that is non-totally-unimodular to file `OUT-MAT`; default: skip computation.

**Advanced options**:
  - `--stats`            Print statistics about the computation to stderr.
  - `--time-limit LIMIT` Allow at most `LIMIT` seconds for the computation.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)
If `IN-MAT` is `-` then the [matrix](\ref file-formats-matrix) is read from stdin.
If `OUT-OPS` or `OUT-MAT` is `-` then the list of operations (resp. the [matrix](\ref file-formats-matrix)) is written to stdout.

### C Interface ###

The corresponding function in the library is

  - CMRctuTest() tests a matrix for being complement totally unimodular.

and is defined in \ref ctu.h.


## Applying Complement Operations ##

The command

    cmr-ctu IN-MAT OUT-MAT [OPTION]...

applies a sequence of row or column complement operations the matrix given in file `IN-MAT` and writes the result to `OUT-MAT`.

**Options**:
  - `-i FORMAT`   Format of file `IN-MAT`; default: [dense](\ref dense-matrix).
  - `-o FORMAT` Format of file `OUT-MAT`; default: same as for `IN-MAT`.
  - `-r ROW`    Apply row complement operation to row `ROW`.
  - `-c COLUMN` Apply column complement operation to column `COLUMN`.

**Advanced options**:
  - `--stats`            Print statistics about the computation to stderr.
  - `--time-limit LIMIT` Allow at most `LIMIT` seconds for the computation.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)
If `IN-MAT` is `-` then the matrix is read from stdin.
If `OUT-MAT` is `-` then the matrix is written to stdout.

## C Interface ##

The corresponding function in the library is

  - CMRctuComplementRowColumn() carries out a row- or column-complement operation for a matrix.

and is defined in \ref ctu.h.

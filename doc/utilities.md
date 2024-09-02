# Basic Utilities # {#utilities}

This tool is useful for basic matrix operations:

  - Transposing the input matrix.
  - Turning a submatrix of a matrix into an explicit matrix.
  - Computing the support matrix of the input matrix.
  - Computing the signed support matrix (negative entries are turned into \f$ -1 \f$'s, positive into \f$ +1 \f$'s) of the input matrix.

## Matrix Utilities ##

The command

    cmr-matrix IN-MAT OUT-MAT [OPTION]...

copies the [matrix](\ref file-formats-matrix) from file `IN-MAT` to file `OUT-MAT`, potentially applying certain operations.

**Options:**
  - `-i FORMAT` Format of file `IN-MAT`; default: [dense](\ref dense-matrix).
  - `-o FORMAT` Format of file `OUT-MAT`; default: same as format of `IN-MAT`.
  - `-S IN-SUB` Consider the [submatrix](\ref file-formats-submatrix) of `IN-MAT` specified in file `IN-SUB` instead of `IN-MAT` itself; can be combined with other operations.
  - `-t`        Transpose the matrix; can be combined with other operations.
  - `-c`        Compute the support matrix instead of copying.
  - `-C`        Compute the signed support matrix instead of copying.
  - `-r`        Randomize the output matrix by randomly permuting rows/columns.
  - `-R2 NUM`   Randomize the output matrix by performing `NUM` binary random pivots.
  - `-R3 NUM`   Randomize the output matrix by performing `NUM` ternary random pivots.
  - `-d`        Use double arithmetic instead of integers.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)
If `IN-MAT` or `IN-SUB` is `-` then the input matrix (resp. submatrix) is read from stdin.
If `OUT-MAT` is `-` then the output matrix is written to stdout.

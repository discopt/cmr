# Matrix Conversions # {#conversions}

This tool is useful for basic matrix operations:

  - Transposing the input matrix.
  - Computing the support matrix of the input matrix.
  - Computing the signed support matrix (negative entries are turned into \f$ -1 \f$'s, positive into \f$ +1 \f$'s) of the input matrix.

## Usage ##

The command

    cmr-convert-matrix IN-MAT OUT-MAT [OPTION]...

copies the matrix from file `IN-MAT` to file `OUT-MAT`, potentially applying certain operations.

**Options:**
  - `-i FORMAT` Format of file `IN-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-o FORMAT` Format of file `OUT-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: same as format of `IN-MAT`.
  - `-t`        Transpose the matrix; can be combined with other operations.
  - `-s`        Compute the support matrix instead of copying.
  - `-S`        Compute the signed support matrix instead of copying.
  - `-d`        Use double arithmetic instead of integers.

If `IN-MAT` is `-` then the input matrix is read from stdin.
If `OUT-MAT` is `-` then the output matrix is written to stdout.

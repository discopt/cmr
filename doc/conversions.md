# Matrix Conversions # {#conversions}

This tool is useful for basic matrix operations:

  - Transposing the input matrix matrix.
  - Computing the support matrix of the input matrix.
  - Computing the signed support matrix (negative entries are turned into \f$ -1 \f$'s, positive into \f$ +1 \f$'s) of the input matrix.
  - Computing a [Camion-signed](\ref camion) version of the input matrix.

## Usage ##

The executable `cmr-convert-matrix` copies the input matrix, potentially applying operations.

    ./cmr-convert-matrix [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output matrices; default: `dense`.
  - `-t` Output the transposed matrix (can be combined with other operations).
  - `-s` Create the support matrix instead of copying.
  - `-s` Create the signed support matrix instead of copying.
  - `-c` Creates the [Camion-signed](\ref camion) version instead of copying.
  - `-d` Use double arithmetic.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.

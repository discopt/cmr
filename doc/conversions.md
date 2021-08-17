# Matrix Conversions # {#conversions}

This tool is useful for basic matrix operations:

  - Transposing the input matrix matrix.
  - Computing the support matrix of the input matrix.
  - Computing a [Camion-signed](\ref camion) version of the input matrix.

## Usage ##

The executable `cmr-convert-matrix` copies the input matrix, potentially applying operations.

    ./cmr-convert-matrix [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output matrices; default: `dense`.
  - `-s` Create the support matrix instead of copying.
  - `-c` Creates the [Camion-signed](\ref camion) version instead of copying.
  - `-t` Output the transposed matrix (can be combined with other operations).
  - `-d` Use double arithmetic.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.

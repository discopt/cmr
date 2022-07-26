# Matrix Entries Inspection # {#entries}

This tool is useful for basic matrix inspection:

  - Test if the input matrix is integer.
  - Test if the input matrix is ternary, i.e, it has only entries in \f$ \{-1,0,+1\} \f$.
  - Test if the input matrix is binary, i.e, it has only entries in \f$ \{0,+1\} \f$.

## Usage ##

The executable `cmr-entries` inspects the matrix.

    ./cmr-entries [OPTION]... MATRIX

Options:
  - `-i FORMAT`  Format of MATRIX file, among {dense, sparse}; default: dense.
  - `-b`         Tests whether the matrix is binary, i.e., has entries in \f$ \{0,+1\} \f$.
  - `-t`         Tests whether the matrix is ternary, i.e., has entries in \f$ \{-1,0,+1\} \f$.
  - `-I`         Tests whether the matrix is integer.
  - `-t EPSILON` Allows rounding of numbers up to tolerance `EPSILON`; default: \f$ 10^{-9} \f$.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If MATRIX is `-`, then the matrix will be read from stdin.

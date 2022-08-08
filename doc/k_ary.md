# Matrix Entries Inspection # {#k_ary}

This tool is useful for basic matrix inspection:

  - Test if the input matrix is integer.
  - Test if the input matrix is ternary, i.e, it has only entries in \f$ \{-1,0,+1\} \f$.
  - Test if the input matrix is binary, i.e, it has only entries in \f$ \{0,+1\} \f$.


## Recognizing Binary and Ternary Matrices ##

The command

    cmr-k-ary IN-MAT [OPTION...]

determines whether the matrix given in file `IN-MAT` is integer (resp. binary or ternary).

**Options**:
  - `-i FORMAT`  Format of file `IN-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-b`         Test whether the matrix is binary, i.e., has entries in \f$ \{0,+1\} \f$.
  - `-t`         Test whether the matrix is ternary, i.e., has entries in \f$ \{-1,0,+1\} \f$.
  - `-I`         Test whether the matrix is integer.
  - `-e EPSILON` Allows rounding of numbers up to tolerance `EPSILON`; default: \f$ 10^{-9} \f$.

If `IN-MAT` is `-` then the matrix is read from stdin.


## Finding Large Binary or Ternary Submatrices ##

The command

    cmr-k-ary IN-MAT -R OUT-SUB [OPTION...]

finds a large binary (resp. ternary) submatrix of the matrix given in file `IN-MAT`.

**Options:**
  - `-i FORMAT`  Format of file `IN-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-b`         Find a large binary submatrix, i.e., one with only entries in \f$ \{0,+1\} \f$.
  - `-t`         Find a large ternary submatrix, i.e., one with only entries in \f$ \{-1,0,+1\} \f$.
  - `-e EPSILON` Allows rounding of numbers up to tolerance `EPSILON`; default: \f$ 10^{-9} \f$.

If `IN-MAT` is `-` then the matrix is read from stdin.
If `OUT-SUB` is `-` then the submatrix is written to stdout.

## Algorithm ##

The implemented algorithm successively removes a row or a column with the maximum number of forbidden entries.

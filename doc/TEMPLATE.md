# Description # {#shortcut}

Definitions, basic tool description.


## Application ##

The command

    cmr-NAME IN-MAT OUT-MAT [OPTION...]

reads the input matrix and outputs a submatrix.
    
**Options**:
  - `-i FORMAT`  Format of file `IN-MAT`; default: [dense](\ref dense-matrix).
  - `-o FORMAT`  Format of file `OUT-MAT`; default: [dense](\ref dense-matrix).
  - `-t`         Consider the transpose of the matrix.
  - `-O OUT-MAT` Write the ... matrix to file `OUT-MAT`; default: stdout.
  - `-N OUT-SUB` Write a minimal .. submatrix to file `OUT-SUB`; default: skip computation.

**Advanced options**:
  - `--stats`            Print statistics about the computation to stderr.
  - `--time-limit LIMIT` Allow at most `LIMIT` seconds for the computation.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)
If `IN-MAT` or `IN-FOO` is `-` then the matrix (resp. foo) is read from stdin.
If `OUT-MAT` is `-` then the matrix is written to stdout.

### Algorithm ###

The implemented algorithm...
For a matrix \f$ M \in \{0,1\}^{m \times n}\f$ with \f$ k \f$ nonzeros it runs in \f$ \mathcal{O}( ? ) \f$ time.

### C Interface ###

The corresponding function in the library is

  - CMRjobOnMatrix() does something.

and is defined in \ref header.h.

...

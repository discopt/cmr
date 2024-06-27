# Series-Parallel Matrices # {#series-parallel}

A matrix \f$ A \in \{-1,0,+1\}^{m \times n} \f$ is called **series-parallel** if it can be obtained from a \f$ 0 \f$-by-\f$ 0 \f$ matrix by successively adjoining

  - a zero row/column vector,
  - a (negated) standard unit row/column vector, or
  - a (negated) copy of an existing row/column.

The removal of such a row/column is called an **SP-reduction**.
A matroid is called **series-parallel** if it is [represented](\ref matroids) by a series-parallel matrix.
This is equivalent to being the graphic matroid of a series-parallel graph.

**Theorem.** A matrix \f$ A \in \{-1,0,1\}^{m \times n} \f$ is *either* series-parallel or it contains, up to scaling of rows/columns with \f$ -1 \f$,

  - a \f$ 2 \f$-by-\f$ 2 \f$ submatrix \f$ M_2 := \begin{pmatrix} -1 & 1 \\ 1 & 1 \end{pmatrix} \f$,
  - a \f$ 3 \f$-by-\f$ 3 \f$ submatrix \f$ M_3' := \begin{pmatrix} 1 & 1 & 0 \\ 1 & 1 & 1 \\ 0 & 1 & 1 \end{pmatrix} \f$, or
  - a \f$ k \f$-by-\f$ k \f$-submatrix \f$ M_k := \begin{pmatrix}
    1 & 0 & 0 & 0 & \dotsb & 0 & 1 \\
    1 & 1 & 0 & 0 & \dotsb & 0 & 0 \\
    0 & 1 & 1 & 0 & \dotsb & 0 & 0 \\
    0 & 0 & 1 & 1 & \ddots & 0 & 0 \\
    \vdots & \vdots & \vdots & \ddots & \ddots & \vdots & \vdots \\
    0 & 0 & 0 & 0 & \dotsb & 1 & 0 \\
    0 & 0 & 0 & 0 & \dotsb & 1 & 1
  \end{pmatrix} \f$ for \f$ k \geq 3 \f$.

The latter two matrices are called **wheel matrices** since their represented matroids are the graphic matroids of wheel graphs.


## Recognizing Series-Parallel Matrices ##

The command

    cmr-series-parallel IN-MAT [OPTION...]

determines whether the [matrix](\ref file-formats-matrix) given in file `IN-MAT` is series-parallel.
If this is not the case, then a maximal number of SP-reductions is carried out, leading to the **reduced** matrix.
Moreover, one can ask for one of the minimal non-series-parallel submatrices above.

**Options:**
  - `-i FORMAT`       Format of file `IN-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-S OUT-SP`       Write the list of series-parallel reductions to file `OUT-SP`; default: skip computation.
  - `-R OUT-REDUCED`  Write the reduced submatrix to file `OUT-REDUCED`; default: skip computation.
  - `-N NON-SUB`      Write a minimal non-series-parallel submatrix to file `NON-SUB`; default: skip computation.
  - `-b`              Test for being binary series-parallel; default: ternary.
  - `-s`              Print statistics about the computation to stderr.

If `IN-MAT` is `-` then the [matrix](\ref file-formats-matrix) is read from stdin.
If `OUT-SP`, `OUT-REDUCED` or `NON-SUB` is `-` then the list of reductions (resp. the [submatrix](\ref file-formats-submatrix)) is written to stdout.

## Algorithm ##

The implemented algorithm is not yet published.
For a matrix \f$ A \in \{0,1\}^{m \times n}\f$ with \f$ k \f$ (sorted) nonzeros it runs in \f$ \mathcal{O}( m + n + k ) \f$ time assuming no hashtable collisions.

## C Interface ##

The corresponding functions in the library are

  - CMRspTestTernary() tests a binary matrix for being series-parallel.
  - CMRspTestBinary() tests a binary matrix for being series-parallel.
  - CMRspDecomposeBinary() tests a binary matrix for being series-parallel, but may also terminate early, returning a 2-separation of \f$ A \f$.
  - CMRspDecomposeTernary() tests a ternary matrix for being series-parallel, but may also terminate early, returning a 2-separation of \f$ A \f$.
  
and are defined in \ref series_parallel.h.

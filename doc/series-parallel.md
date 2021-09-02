# Series-Parallel Matroids # {#series-parallel}

A matrix \f$ A \in \{-1,0,+1\}^{m \times n} \f$ is called **series-parallel** if it can be obtained from a \f$ 0 \f$-by-\f$ 0 \f$ matrix by successively adjoining

  - a zero row/column vector,
  - a (negated) standard unit row/column vector, or
  - a (negated) copy of an existing row/column.

The removal of such a row/column is called an **SP-reduction**.
A matroid is called **series-parallel** if it is represented by a series-parallel matrix.
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

The latter two matrices are called **wheel matrices** since they represent wheel graphs.

## Usage ##

The executable `cmr-series-parallel` tests whether a given matrix \f$ A \f$ is series-parallel.
If this is not the case, then a maximal number of SP-reductions is carried out, leading to the **reduced** matrix.
Moreover, one can ask for one of the minimal non-series-parallel submatrices above.

    ./cmr-series-parallel [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output matrices; default: `dense`.
  - `-sp`       Output the list of series-parallel reductions.
  - `-r`        Output the elements of the reduced matrix.
  - `-R`        Output the reduced matrix.
  - `-n`        Output the elements of a minimal non-series-parallel submatrix.
  - `-N`        Output a minimal non-series-parallel submatrix.

Formats for matrices are \ref dense-matrix, \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.

## Algorithm ##

The implemented algorithm is not yet published.
For a matrix \f$ A \in \{0,1\}^{m \times n}\f$ with \f$ k \f$ (sorted) nonzeros it runs in \f$ \mathcal{O}( m + n + k ) \f$ time assuming no hashtable collisions.

## C Interface ##

  - CMRtestTernarySeriesParallel() tests a binary matrix for being series-parallel.
  - CMRtestBinarySeriesParallel() tests a binary matrix for being series-parallel.
  - CMRdecomposeBinarySeriesParallel() tests a binary matrix for being series-parallel, but may also terminate early, returning a 2-separation of \f$ A \f$.
  - CMRdecomposeTernarySeriesParallel() tests a ternary matrix for being series-parallel, but may also terminate early, returning a 2-separation of \f$ A \f$.

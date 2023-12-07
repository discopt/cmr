# Totally Unimodular Matrices # {#tu}

A matrix \f$ M \in \mathbb{Z}^{m \times n} \f$ is **totally unimodular** if all its square submatrices have a determinant in \f$ \{-1,0,+1\} \f$.
Here, a submatrix does not need to be contiguous, i.e., the matrix \f$M = \begin{pmatrix} 1 & 0 & -1 \\ 1 & 0 & 1 \end{pmatrix} \f$ is not totally unimodular since the submatrix indexed by rows \f$ \{1, 2 \} \f$ and columns \f$ \{ 1,3 \} \f$ has determinant 2.
In particular, every totally unimodular matrix has only entries in \f$ \{-1,0,+1\} \f$ as these are the 1-by-1 submatrices.


## Recognizing Totally Unimodular Matrices  ##

The command

    cmr-tu IN-MAT [OPTION...]

determines whether the matrix given in file `IN-MAT` is totally unimodular.

**Options:**
  - `-i FORMAT`  Format of file `IN-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-D OUT-DEC` Write a decomposition tree of the underlying regular matroid to file `OUT-DEC`; default: skip computation.
  - `-N NON-SUB` Write a minimal non-totally-unimodular submatrix to file `NON-SUB`; default: skip computation.
  - `-s`         Print statistics about the computation to stderr.

**Advanced options:**
  - `--no-direct-graphic`  Check only 3-connected matrices for regularity.
  - `--no-series-parallel` Do not allow series-parallel operations in decomposition tree.

If `IN-MAT` is `-` then the matrix is read from stdin.
If `OUT-DEC` or `NON-SUB` is `-` then the decomposition tree (resp. the submatrix) is written to stdout.

## Algorithms ##

The implemented default recognition algorithm is based on [Implementation of a unimodularity test](https://doi.org/10.1007/s12532-012-0048-x) by Matthias Walter and Klaus Truemper (Mathematical Programming Computation, 2013).
It first runs \ref camion to reduce the question to that of [recognizing regular matroids](\ref regular).
Please cite the paper in case the implementation contributed to your research:

    @Article{WalterT13,
      author    = {Walter, Matthias and Truemper, Klaus},
      title     = {Implementation of a unimodularity test},
      journal   = {Mathematical Programming Computation},
      year      = {2013},
      volume    = {5},
      number    = {1},
      pages     = {57--73},
      issn      = {1867-2949},
      doi       = {10.1007/s12532-012-0048-x},
      publisher = {Springer-Verlag},
    }

In order to repeat experiments described in the paper above, the function can be parameterized as to use exponential-time algorithms.
The first is based on the criterion of Ghouila-Houri and runs in time \f$ \mathcal{O}( (m + n) 3^{\min(m, n)}) \f$.
The second enumerates square [Eulerian submatrices](https://www.ams.org/journals/proc/1965-016-05/S0002-9939-1965-0180568-2/) and runs in time \f$ \mathcal{O}( (m+n) 2^{ m + n } ) \f$.

## C Interface ##

The corresponding function in the library is

  - CMRtestTotalUnimodularity() tests a matrix for being totally unimodular.

and is defined in \ref tu.h.

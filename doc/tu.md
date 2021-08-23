# Totally Unimodular Matrices # {#tu}

A matrix \f$ M \in \mathbb{Z}^{m \times n} \f$ is **totally unimodular** if all its square submatrices have a determinant in \f$ \{-1,0,+1\} \f$.
Here, a submatrix does not need to be contiguous, i.e., the matrix \f$M = \begin{pmatrix} 1 & 0 & -1 \\ 1 & 0 & 1 \end{pmatrix} \f$ is not totally unimodular since the submatrix indexed by rows \f$ \{1, 2 \} \f$ and columns \f$ \{ 1,3 \} \f$ has determinant 2.
In particular, every totally unimodular matrix has only entries in \f$ \{-1,0,+1\} \f$ as these are the 1-by-1 submatrices.

## Usage ##

The executable `cmr-tu` determines whether a given matrix \f$ M \f$ is totally unimodular.

    ./cmr-tu [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output matrices; default: `dense`.
  - `-d`        Output the decomposition tree of the underlying regular matroid.
  - `-s`        Output the elements of a minimal non-totally-unimodular submatrix.
  - `-S`        Output a minimal non-totally-unimodular submatrix.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.


## Algorithm ##

The implemented recognition algorithm is based on [Implementation of a unimodularity test](https://doi.org/10.1007/s12532-012-0048-x) by Matthias Walter and Klaus Truemper (Mathematical Programming Computation, 2013).
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

## C Interface ##

The functionality is defined in \ref tu.h.
The main functions are:

  - CMRtestTotalUnimodularity() tests a matrix for being totally unimodular.

# Graphic / Cographic / Planar Matroids# {#graphic}

Let \f$ G = (V,E) \f$ be a graph and let \f$ T \f$ be a spanning forest.
The matrix \f$ M(G,T) \in \{0,1\}^{T \times (E \setminus T)} \f$ defined via
\f[
  M(G,T)_{e,f} := \begin{cases}
    1 & \text{if } e \text{ is contained in the unique cycle of } T \cup \{f\} \\
    0 & \text{otherwise}
  \end{cases}
\f]
is called the **graphic matrix** of \f$ G \f$ with respect to \f$ T \f$.
A binary matrix \f$ M \f$ is called **graphic** if there exists a graph \f$ G \f$ with a spanning forest \f$ T \f$ such that \f$ M = M(G,T) \f$.
Moreover, \f$ M \f$ is called **cographic** if \f$ M^{\textsf{T}} \f$ is graphic, and 
it is called **planar** if it is graphic and cographic.

## Usage ##

The executable `cmr-graphic` converts graphs to graphic matrices and vice versa.
In particular, for a given matrix \f$ M \f$, it determines whether \f$ M \f$ is (co)graphic.

    ./cmr-graphic [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output; default: `edgelist` if input is a matrix and `dense` if input is a graph.
  - `-t`        Tests for being / converts to cographic matrix.

Formats for matrices are \ref dense-matrix, \ref sparse-matrix.
Formats for graphs are \ref edge-list and \ref dot.
If FILE is `-`, then the input will be read from stdin.

## Algorithm ##

The implemented recognition algorithm is based on [An Almost Linear-Time Algorithm for Graph Realization](https://doi.org/10.1287/moor.13.1.99) by Robert E. Bixby and Donald K. Wagner (Mathematics of Operations Research, 1988).
For a matrix \f$ M \in \{0,1\}^{m \times n}\f$ with \f$ k \f$ nonzeros it runs in \f$ \mathcal{O}( k \cdot \alpha(k, m) ) \f$ time, where \f$ \alpha(\cdot) \f$ denotes the inverse Ackerman function.


## C Interface ##

The functionality is defined in \ref graphic.h.
The main functions are:

  - CMRcomputeGraphicMatrix() constructs a graphic matrix for a given graph.
  - CMRtestGraphicMatrix() tests a matrix for being graphic.
  - CMRtestCographicMatrix() tests a matrix for being cographic.

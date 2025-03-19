# Graphic / Cographic / Planar Matrices # {#graphic}

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


## Recognizing Graphic Matrices ##

The command

    cmr-graphic IN-MAT [OPTION]...

determines whether the [matrix](\ref file-formats-matrix) given in file `IN-MAT` is (co)graphic.

**Options**:
  - `-i FORMAT`    Format of file `IN-MAT`; default: [dense](\ref dense-matrix).
  - `-t`           Test for being cographic; default: test for being graphic.
  - `-G OUT-GRAPH` Write a graph to file `OUT-GRAPH`; default: skip computation.
  - `-T OUT-TREE`  Write a spanning tree to file `OUT-TREE`; default: skip computation.
  - `-D OUT-DOT`   Write a dot file `OUT-DOT` with the graph and the spanning tree; default: skip computation.

**Advanced options**:
  - `--stats`            Print statistics about the computation to stderr.
  - `--time-limit LIMIT` Allow at most `LIMIT` seconds for the computation.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

If `IN-MAT` is `-` then the [matrix](\ref file-formats-matrix) is read from stdin.

If `OUT-GRAPH`, `OUT-TREE`, `OUT-DOT` or `NON-SUB` is `-` then the graph (resp. the tree, dot file or non-(co)graphic [submatrix](\ref file-formats-submatrix)) is written to stdout.

### Algorithm ###

The implemented recognition algorithm is based on [An Almost Linear-Time Algorithm for Graph Realization](https://doi.org/10.1287/moor.13.1.99) by Robert E. Bixby and Donald K. Wagner (Mathematics of Operations Research, 1988).
For a matrix \f$ M \in \{0,1\}^{m \times n}\f$ with \f$ k \f$ nonzeros it runs in \f$ \mathcal{O}( k \cdot \alpha(k, m) ) \f$ time, where \f$ \alpha(\cdot) \f$ denotes the inverse Ackerman function.

### C Interface ###

The corresponding functions in the library are

  - CMRgraphicTestMatrix() tests a matrix for being graphic.
  - CMRgraphicTestTranspose() tests a matrix for being cographic.

and are defined in \ref network.h.


## Computing graphic matrices ##

The command

    cmr-graphic -c IN-GRAPH OUT-MAT [OPTION]...

computes a (co)graphic [matrix](\ref file-formats-matrix) corresponding to the graph from file `IN-GRAPH` and writes it to `OUT-MAT`.

**Options**:
  - `-o FORMAT`    Format of file `OUT-MAT`; default: [dense](\ref dense-matrix).
  - `-t`           Return the transpose of the graphic matrix.
  - `-T IN-TREE`   Read a tree from file `IN-TREE`; default: use first specified arcs as tree edges.

**Advanced options**:
  - `--stats`            Print statistics about the computation to stderr.
  - `--time-limit LIMIT` Allow at most `LIMIT` seconds for the computation.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

If `IN-GRAPH` or `IN-TREE` is `-` then the graph (resp. tree) is read from stdin.

If `OUT-MAT` is `-` then the [matrix](\ref file-formats-matrix) is written to stdout.

### C Interface ###

The corresponding function in the library is

  - CMRgraphicComputeMatrix() constructs a graphic matrix for a given graph.

and is defined in \ref network.h.

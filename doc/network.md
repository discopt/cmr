# Network Matrices # {#network}

Let \f$ D = (V,A) \f$ be a digraph and let \f$ T \f$ be an (arbitrarily) directed spanning forest of the underlying undirected graph.
The matrix \f$ M(D,T) \in \{-1,0,1\}^{T \times (A \setminus T)} \f$ defined via
\f[
  M(D,T)_{a,(v,w)} := \begin{cases}
    +1 & \text{if the unique $v$-$w$-path in $T$ passes through $a$ forwardly}, \\
    -1 & \text{if the unique $v$-$w$-path in $T$ passes through $a$ backwardly}, \\
    0  & \text{otherwise}
  \end{cases}
\f]
is called the **network matrix** of \f$ D \f$ with respect to \f$ T \f$.
A matrix \f$ M \f$ is called **network matrix** if there exists a digraph \f$ D \f$ with a directed spanning forest \f$ T \f$ such that \f$ M = M(D,T) \f$.
Moreover, \f$ M \f$ is called **conetwork matrix** if \f$ M^{\textsf{T}} \f$ is a network matrix.


## Recognizing Network Matrices ##

The command

    cmr-network IN-MAT [OPTION]...

determines whether the matrix given in file `IN-MAT` is (co)network.

**Options**:
  - `-i FORMAT`    Format of file `IN-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-t`           Test for being conetwork; default: test for being network.
  - `-G OUT-GRAPH` Write a digraph to file `OUT-GRAPH`; default: skip computation.
  - `-T OUT-TREE`  Write a directed spanning tree to file `OUT-TREE`; default: skip computation.
  - `-D OUT-DOT`   Write a dot file `OUT-DOT` with the digraph and the directed spanning tree; default: skip computation.
  - `-N NON-SUB`   Write a minimal non-network submatrix to file `NON-SUB`; default: skip computation.
  - `-s`           Print statistics about the computation to stderr.

If `IN-MAT` is `-` then the matrix is read from stdin.
If `OUT-GRAPH`, `OUT-TREE`, `OUT-DOT` or `NON-SUB` is `-` then the graph (resp. the tree, dot file or non-(co)network submatrix) is written to stdout.

### Algorithm ###

The implemented recognition algorithm first tests the support matrix of \f$ M \f$ for being [(co)graphic](\ref graphic) and uses \ref camion for testing whether \f$ M \f$ is signed correctly.

### C Interface ###

The corresponding functions in the library are

  - CMRtestNetworkMatrix() tests a matrix for being network.
  - CMRtestConetworkMatrix() tests a matrix for being conetwork.

and are defined in \ref network.h.


## Computing network matrices ##

The command

    cmr-network -c IN-GRAPH OUT-MAT [OPTION]...

computes a (co)network matrix corresponding to the digraph from file `IN-GRAPH` and writes it to `OUT-MAT`.

**Options**:
  - `-o FORMAT`    Format of file `OUT-MAT`, among `dense` for \ref dense-matrix and `sparse` for \ref sparse-matrix; default: dense.
  - `-t`           Return the transpose of the network matrix.
  - `-T IN-TREE`   Read a directed tree from file `IN-TREE`; default: use first specified arcs as tree edges.
  - `-s`           Print statistics about the computation to stderr.

If `IN-GRAPH` or `IN-TREE` is `-` then the digraph (resp. directed tree) is read from stdin.
If `OUT-MAT` is `-` then the matrix is written to stdout.

### C Interface ###

The corresponding function in the library is

  - CMRcomputeNetworkMatrix() constructs a network matrix for a given digraph.

and is defined in \ref network.h.

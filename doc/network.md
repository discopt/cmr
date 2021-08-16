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
A matrix \f$ A \f$ is called **network matrix** if there exists a digraph \f$ D \f$ with a directed spanning forest \f$ T \f$ such that \f$ A = M(D,T) \f$.
Moreover, \f$ A \f$ is called **conetwork matrix** if \f$ A^{\textsf{T}} \f$ is a network matrix.

## Usage ##

The executable `cmr-network` converts digraphs to (co)network matrices and vice versa.
In particular, for a given matrix \f$ A \f$, it determines whether \f$ A \f$ is (co)network.

    ./cmr-network [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output; default: `edgelist` if input is a matrix and `dense` if input is a graph.
  - `-t`        Tests for being / converts to conetwork matrix.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
Formats for graphs are \ref edge-list and \ref dot (output-only).
If FILE is `-`, then the input will be read from stdin.

## Algorithm ##

The implemented algorithm is based on [An Almost Linear-Time Algorithm for Graph Realization](https://doi.org/10.1287/moor.13.1.99) by Robert E. Bixby and Donald K. Wagner (Mathematics of Operations Research, 1988).
For a matrix \f$ A \in \{0,1\}^{m \times n}\f$ with \f$ k \f$ nonzeros it runs in \f$ \mathcal{O}( k \cdot \alpha(k, m) ) \f$ time, where \f$ \alpha(\cdot) \f$ denotes the inverse Ackerman function.

## C Interface ##

The functionality is defined in \ref network.h.
The main functions are:

  - CMRcomputeNetworkMatrix() constructs a network matrix for a given digraph.
  - CMRcomputeConetworkMatrix() constructs a conetwork matrix for a given digraph.
  - CMRtestNetworkMatrix() tests a matrix for being network.
  - CMRtestConetworkMatrix() tests a matrix for being conetwork.

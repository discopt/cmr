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

## Usage ##

The executable `cmr-network` converts digraphs to (co)network matrices and vice versa.
In particular, for a given matrix \f$ M \f$, it determines whether \f$ M \f$ is (co)network.

    ./cmr-network [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output; default: `edgelist` if input is a matrix and `dense` if input is a graph.
  - `-t`        Tests for being / converts to conetwork matrix.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
Formats for graphs are \ref edge-list and \ref dot (output-only).
If FILE is `-`, then the input will be read from stdin.

## Algorithm ##

The implemented recognition algorithm first tests the support matrix of \f$ M \f$ for being [(co)graphic](\ref graphic) and uses \ref camion for testing whether \f$ M \f$ is signed correctly.

## C Interface ##

The functionality is defined in \ref network.h.
The main functions are:

  - CMRcomputeNetworkMatrix() constructs a network matrix for a given digraph.
  - CMRtestNetworkMatrix() tests a matrix for being network.
  - CMRtestConetworkMatrix() tests a matrix for being conetwork.

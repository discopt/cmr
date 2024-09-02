# Instance Generators # {#generators}

The library comes with several instance generators.
However, they are disabled for compilation by default, and can be enabled by adding

    -DGENERATORS=on

to the cmake call.
In the following we provide details for each of these generators.

## Graphic Matrices ##

The executable `cmr-generate-graphic` creates a random \f$ m \f$-by-\f$n\f$ matrix that is graphic and writes it to `stdout`.
To this end, it first constructs a spanning tree with \f$ m+1 \f$ nodes by subsequently connecting each new leaf to one of the previously created nodes.
The new edge then indexes a row of the matrix.
Second, \f$ n \f$ random pairs of distinct nodes are generated that are then connected by a new edge indexed by a column.
The matrix is the representation matrix of the final graph with respect to the spanning tree.
It can be called as follows, where \f$ m \f$ equals ROWS and \f$ n \f$ equals COLS.

    ./cmr-generate-graphic ROWS COLS [OPTION]...

**Options:**
  - `-B NUM`    Benchmarks the recognition algorithm for the created matrix with NUM repetitions.
  - `-o FORMAT` Format of output; default: [dense](\ref dense-matrix).

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

## Network Matrices ##

The executable `cmr-generate-network` creates a random \f$ m \f$-by-\f$n\f$ network matrix, either ternary (default) or binary and writes it to `stdout`.
In order to generate a ternary matrix, it first constructs a random graphic matrix (see above) and then modifies its signs via \ref camion.
For the binary matrix, a random spanning tree is constructed and oriented to become an arborescence with root.
For the co-tree edges, only node pairs are considered, where the second node lies on the directed path from the first node to the root.
It can be called as follows, where \f$ m \f$ equals ROWS and \f$ n \f$ equals COLS.

    ./cmr-generate-network ROWS COLS [OPTION]...

**Options:**
  - `-b`        Restrict to binary network matrices based on arborescences.
  - `-B NUM`    Benchmarks the recognition algorithm for the created matrix with NUM repetitions.
  - `-o FORMAT` Format of output; default: [dense](\ref dense-matrix).

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

## Random Matrices ##

The executable `cmr-generate-random` creates a random \f$ m \f$-by-\f$n\f$ binary matrix whose entries are chosen uniformly at random with a given probability \f$ p \f$ and writes it to `stdout`.
It can be called as follows, where \f$ m \f$ equals `ROWS` and \f$ n \f$ equals `COLS`.

    ./cmr-generate-random ROWS COLS p [OPTION]...

**Options:**
  - `-o FORMAT` Format of output; default: [dense](\ref dense-matrix).

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

## Random Perturbations ##

The executable `cmr-perturb-random` copies the matrix from file `IN-MAT` to the file `OUT-MAT` after applying random perturbations.
It can be called as follows.

    ./cmr-perturb-random IN-MAT OUT-MAT [OPTION]...

**Options:**
  - `-i FORMAT` Format of file `IN-MAT`; default: dense.
  - `-o FORMAT` Format of file `OUT-MAT`; default: same as format of `IN-MAT`.
  - `-0 NUM`    Turn NUM randomly chosen nonzero entries to 0s.
  - `-1 NUM`    Turn NUM randomly chosen zero entries into 1s.
  - `--1 NUM`   Turn NUM randomly chosen zero entries into -1s.
  - `-b NUM`    Flip NUM randomly chosen entries over the binary field.
  - `-t NUM`    Flip NUM randomly chosen entries over the ternary field.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)
If `IN-MAT` is `-`, then the input matrix is read from stdin.
If `OUT-MAT` is `-`, then the output matrix is written to stdout.

## Wheel Matrices ##

The executable `cmr-generate-wheel` creates an \f$ n \f$-by-\f$ n \f$ **wheel matrix** \f$ W_n \f$ or a variant \f$ W'_n \f$ and writes it to `stdout`.

\f[
  W_n = 
  \begin{pmatrix}
    1 & 1 & 0 & 0 & 0 & \dots & 0 & 0 \\
    0 & 1 & 1 & 0 & 0 & \dots & 0 & 0 \\
    0 & 0 & 1 & 1 & 0 & \dots & 0 & 0 \\
    \vdots & & & \ddots & \ddots & & & \vdots \\
    \vdots & & & & \ddots & \ddots & & \vdots \\
    0 & 0 & \dots & 0 & 0 & 1 & 1 & 0 \\
    0 & 0 & \dots & 0 & 0 & 0 & 1 & 1 \\
    1 & 0 & \dots & 0 & 0 & 0 & 0 & 1
  \end{pmatrix},
  \qquad
  \qquad
  \qquad
  W'_n = 
  \begin{pmatrix}
    1 & 1 & 0 & \color{red}{1} & \color{red}{1} & \dots & \color{red}{1} & \color{red}{1} \\
    0 & 1 & 1 & \color{red}{1} & \color{red}{1} & \dots & \color{red}{1} & \color{red}{1} \\
    0 & 0 & 1 & 1 & 0 & \dots & 0 & 0 \\
    \vdots & & & \ddots & \ddots & & & \vdots \\
    \vdots & & & & \ddots & \ddots & & \vdots \\
    0 & 0 & \dots & 0 & 0 & 1 & 1 & 0 \\
    0 & 0 & \dots & 0 & 0 & 0 & 1 & 1 \\
    1 & 0 & \dots & 0 & 0 & 0 & 0 & 1
  \end{pmatrix}
\f]

It can be called as follows, where \f$ n \f$ equals `ORDER`.

    ./cmr-generate-wheel ORDER [OPTION]...

**Options:**
  - `-01`       In each column with two `0`s in rows 1 and 2, replace them by `1`s; default: off
  - `-o FORMAT` Format of output; default: [dense](\ref dense-matrix).

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

\note
The matrices \f$ W_n \f$ are [graphic](\ref graphic) for all \f$ n \f$.
The modified matrices \f$ W'_n \f$ for odd \f$ n \f$ are **minimally non-[totally unimodular](\ref tu)**, which means that the matrix itself has determinant \f$ |\det W'_n| = 2 \f$, but all of its proper submatrices are totally unimodular.

## Gurobi Coefficient Matrix ##

The executable `cmr-extract-gurobi` extracts the coefficient matrix of a mixed-integer program file that can be read by the [Gurobi solver](https://www.gurobi.com) and writes it to `stdout`.
Building it requires that Gurobi is found by cmake.
To this end, you may need to set the cmake option `-DGUROBI_DIR` to your Gurobi installation directory.
It can be called as follows.

    ./cmr-extract-gurobi MIPFILE [OPTION]...

Options:
  - `-o FORMAT` Format of output matrix; default: `dense`.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)
`MIPFILE` must refer to a file that Gurobi can read.


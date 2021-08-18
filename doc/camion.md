# Camion's Signing Algorithm # {#camion}

A key tool for the recognition of [network matrices](\ref network), [totally unimodular matrices](\ref tu) and [balanceable matrices](\ref balanced) is Camion's signing algorithm.
Its input is a binary matrix \f$ M \in \{0,1\}^{m \times n} \f$ and it computes a ternary matrix \f$ M' \in \{-1,0,+1\}^{m \times n} \f$ with the same support, i.e., the same nonzero entries.
The output matrix \f$ M' \f$ is guaranteed to be [balanced](\ref balanced) if (and only if) \f$ M \f$ was [balanceable](\ref balanced).
This means that for every square submatrix of \f$ M' \f$ with two nonzero entries per row and per column the the sum of all entries is divisible by 4.
The algorithm always outputs a **Camion-signed matrix**, regardless of whether the input matrix was balanceable.
In particular, it does not recognize whether it deals with a balanceable matrix or not.

If \f$ M \f$ is balanceable, then \f$ M' \f$ is unique up to scaling rows/columns with \f$ -1 \f$.


## Usage ##

The executable `cmr-camion` checks whether a given matrix \f$ M \f$ is Camion-signed.

    ./cmr-camion [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-o FORMAT` Format of output matrices; default: `dense`.
  - `-s`        Output the elements of a minimal non-camion submatrix.
  - `-S`        Output a minimal non-camion submatrix.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.

\note A [Camion-signed](\ref camion) version of a matrix can be **computed** via [cmr-convert-matrix -c](\ref conversions).

## Algorithm ##

The implemented recognition algorithm is based on Section 18 of [A decomposition theory for matroids. V. Testing of matrix total unimodularity](https://doi.org/10.1016/0095-8956(90)90030-4) by Klaus Truemper (Journal of Combinatorial Theory, Series B, 1990).
For a matrix \f$ M \in \{-1,0,1\}^{m \times n} \f$ it runs in \f$ \mathcal{O}( \min(m^2 \cdot n, m \cdot n^2) ) \f$ time.

## C Interface ##

The functionality is defined in \ref camion.h.
The main functions are:

  - CMRtestCamionSigned() tests a matrix for being Camion-signed.

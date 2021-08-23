# (Strongly) k-Modular and Unimodular Matrices # {#k-modular}

Consider a matrix \f$ M \in \mathbb{Z}^{m \times n} \f$ of rank \f$ r \f$.
The matrix \f$ M \f$ is called **k-modular** (for some \f$ k \in \{1,2,\dotsc \} \f$) if these two conditions are satisfied:

  - for some column basis \f$ B \subseteq \{1,2,\dotsc,n\} \f$ of \f$ M \f$, the greatest common divisor of the determinants of all \f$r \f$-by-\f$ r \f$ submatrices of \f$ M_{\star,B} \f$ is equal to \f$ k \f$.
  - The matrix \f$ X \f$ such that \f$ M = M_{\star,B} X \f$ is [totally unimodular](\ref tu).

In case \f$ M \f$ has full row-rank, the first property requires that the determinant of any basis matrix shall be \f$ k \f$, while the second property requires that \f$ M_{\star,B}^{-1} M \f$ is [totally unimodular](\ref tu).
Otherwise, \f$ M_{\star,B} \f$ is singular, and hence the property is more technical.

\note k-modularity is independent of the choice of the column basis \f$ B \f$.

Additionally, \f$ M \f$ is called **strongly k-modular** if \f$ M \f$ and \f$ M^{\textsf{T}} \f$ are k-modular.
The special cases with \f$ k = 1 \f$ is called **unimodular** and **strongly unimodular**, respectively.

## Usage ##

The executable `cmr-k-modular` determines whether a given matrix \f$ M \f$ is \f$k\f$-modular (and determines \f$ k \f$).

    ./cmr-k-modular [OPTION]... FILE

Options:
  - `-i FORMAT` Format of input FILE; default: `dense`.
  - `-t`        Test \f$ M^{\textsf{T}} \f$ instead.
  - `-s`        Test for strong \f$ k \f$-modularity, i.e., test \f$ M \f$ and \f$ M^{\textsf{T}} \f$.
  - `-u`        Test only for unimodularity, i.e., \f$ 1 \f$-modularity.

Formats for matrices are \ref dense-matrix and \ref sparse-matrix.
If FILE is `-`, then the input will be read from stdin.

## C Interface ##

The functionality is defined in \ref k_modular.h.
The main functions are:

  - CMRtestKmodularity() tests a matrix for being \f$k\f$-modular.
  - CMRtestStrongKmodularity() tests a matrix for being strongly \f$k\f$-modular.
  - CMRtestUnimodularity() tests a matrix for being unimodular.
  - CMRtestStrongUnimodularity() tests a matrix for being strongly unimodular.

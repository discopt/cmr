# (Strongly) Equimodular and Unimodular Matrices # {#equimodular}

Consider a matrix \f$ M \in \mathbb{Z}^{m \times n} \f$ of rank \f$ r \f$.
The matrix \f$ M \f$ is called **equimodular** *with determinant gcd* \f$ k \in \{1,2,\dotsc \} \f$ if these two conditions are satisfied:

  - for some column basis \f$ B \subseteq \{1,2,\dotsc,n\} \f$ of \f$ M \f$, the greatest common divisor of the determinants of all \f$r \f$-by-\f$ r \f$ submatrices of \f$ M_{\star,B} \f$ is equal to \f$ k \f$.
  - The matrix \f$ X \f$ such that \f$ M = M_{\star,B} X \f$ is [totally unimodular](\ref tu).

In case \f$ M \f$ has full row-rank, the first property requires that the determinant of any basis matrix shall be \f$ \pm k \f$, while the second property requires that \f$ M_{\star,B}^{-1} M \f$ is [totally unimodular](\ref tu).
Otherwise, \f$ M_{\star,B} \f$ is not square, and hence the property is more technical.

\note Equimodularity is independent of the choice of the column basis \f$ B \f$.

Additionally, \f$ M \f$ is called **strongly equimodular** if \f$ M \f$ and \f$ M^{\textsf{T}} \f$ are both equimodular, which implies that they are equimodular for the same gcd determinants.
The special cases with \f$ k = 1 \f$ are called **unimodular** and **strongly unimodular**, respectively.

## Usage ##

The executable `cmr-equimodular` determines whether a given [matrix](\ref file-formats-matrix) \f$ M \f$ with determinant gcd \f$ k \f$.

    ./cmr-equimodular IN-MAT [OPTION]...

**Options:**
  - `-i FORMAT` Format of of file `IN-MAT`; default: [dense](\ref dense-matrix).
  - `-t`        Test \f$ M^{\textsf{T}} \f$ instead.
  - `-s`        Test for strong equimodularity.
  - `-u`        Test only for unimodularity, i.e., \f$ k = 1 \f$.

**Advanced options**:
  - `--stats`            Print statistics about the computation to stderr.
  - `--time-limit LIMIT` Allow at most `LIMIT` seconds for the computation.

Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

If `IN-MAT` is `-`, then the input will be read from stdin.

## C Interface ##

The functionality is defined in \ref equimodular.h.
The main functions are:

  - CMRequimodularTest() tests a matrix for being equimodular.
  - CMRequimodularTestStrong() tests a matrix for being strongly equimodular.
  - CMRunimodularTest() tests a matrix for being unimodular.
  - CMRunimodularTestStrong() tests a matrix for being strongly unimodular.

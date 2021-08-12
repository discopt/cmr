# Totally Unimodular Matrices # {#tu}

A matrix \f$ A \in \mathbb{Z}^{m \times n} \f$ is **totally unimodular** if all its square submatrices have a determinant in \f$ \{-1,0,+1\} \f$.
Here, a submatrix does not need to be contiguous, i.e., the matrix \f$A = \begin{pmatrix} 1 & 0 & -1 \\ 1 & 0 & 1 \end{pmatrix} \f$ is not totally unimodular since the submatrix indexed by rows \f$ \{1, 2 \} \f$ and columns \f$ \{ 1,3 \} \f$ has determinant 2.
In particular, \f$ A \in \{-1,0,+1\} \f$ holds due the 1-by-1 submatrices.

## Usage ##



# Representation of Matroids # {#matroids}

The matroid **represented** by a matrix \f$ M \in \mathbb{F}^{m \times n} \f$ has \f$ m + n \f$ elements that correspond to the columns of \f$ [ \mathbb{I} \mid M ] \f$.
A subset of these columns is independent if and only if it is linearly independent over the field \f$ \mathbb{F} \f$.
Applying a **pivot operation** (over \f$ \mathbb{F} \f$) and exchanging the corresponding row and column elements leaves
the matroid unchanged.
In CMR, pivots are carried out via \ref CMRchrmatBinaryPivot, \ref CMRchrmatBinaryPivots, \ref CMRchrmatTernaryPivot and \ref CMRchrmatTernaryPivots.

## Minors ##

A **minor** of a represented matroid is another one obtained by **deleting** or **contracting** elements.
An element associated with a row (resp. column) can be contracted (resp. deleted) by removing the corresponding matrix row (resp. column).
In order to delete a row element or contract a column element one must first pivot such that the element type (row/column) is changed.

A minor of a matroid represented by matrix \f$ M \f$ is represented by means of a \ref CMR_MINOR object.
It consists of a (potentially empty) array of pivots and a \ref CMR_SUBMAT object indicating a submatrix \f$ M' \f$ of the matrix obtained from \f$ M \f$ after applying the pivots.
Moreover, the **type** field indicates a certain structure of \f$ M' \f$:

 - A [determinant](\ref CMR_MINOR_TYPE_DETERMINANT) type indicates that \f$ M' \f$ is a submatrix of \f$ M \f$ with \f$ |\det(M')| \geq 2 \f$.
   In particular, no pivots are applied.
 - A [Fano](\ref CMR_MINOR_TYPE_FANO) type indicates that \f$ M' \f$ represents the Fano matroid \f$ F_7 \f$.
 - A [Fano-dual](\ref CMR_MINOR_TYPE_FANO_DUAL) type indicates that \f$ M' \f$ represents the dual \f$ F_7^\star \f$ of the Fano matroid.
 - A [K5](\ref CMR_MINOR_TYPE_K5) type indicates that \f$ M' \f$ represents the graphic matroid \f$ M(K_5) \f$ of the complete graph \f$ K_5 \f$.
 - A [K5-dual](\ref CMR_MINOR_TYPE_K5_DUAL) type indicates that \f$ M' \f$ represents the dual matroid \f$ M(K_5)^\star \f$ of the complete graph \f$ K_5 \f$.
 - A [K33](\ref CMR_MINOR_TYPE_K33) type indicates that \f$ M' \f$ represents the graphic matroid \f$ M(K_{3,3}) \f$ of the complete bipartite graph \f$ K_{3,3} \f$.
 - A [K33-dual](\ref CMR_MINOR_TYPE_K33_DUAL) type indicates that \f$ M' \f$ represents the dual matroid \f$ M(K_{3,3})^\star \f$ of the complete bipartite graph \f$ K_{3,3} \f$.

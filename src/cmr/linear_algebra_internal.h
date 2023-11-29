#ifndef CMR_LINALG_INTERNAL_H
#define CMR_LINALG_INTERNAL_H

#include <cmr/linear_algebra.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Transforms \p matrix into a new matrix by applying integer row operations and row- and column swaps to find
 *        an upper-diagonal basis matrix.
 *
 * The rank \f$ r \f$ is stored in \p *prank and the row- and column permutations are stored in \p *ppermutations,
 * such that the first \f$ r \f$ rows and columns of the resulting matrix form an invertible upper-diagonal matrix.
 * If \p invert is \c true then in this \f$ r \f$-by-\f$ r \f$ submatrix, the largest (in terms of absolute value)
 * entry in each column is on the diagonal.
 */

CMR_ERROR CMRintmatComputeUpperDiagonal(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_INTMAT* matrix,         /**< A matrix */
  bool invert,                /**< Whether the transformed basis columns shall be strictly diagonally dominant. */
  size_t* prank,              /**< Pointer for storing the rank of the basis matrix. */
  CMR_SUBMAT** ppermutations, /**< Pointer for storing the row- and column permutations applied to \p matrix. */
  CMR_INTMAT** presult,       /**< Pointer for storing the resulting int matrix (may be \c NULL). */
  CMR_INTMAT** ptranspose     /**< Pointer for storing the transpose of the result int matrix (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_LINALG_INTERNAL_H */

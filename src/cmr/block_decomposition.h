#ifndef CMR_BLOCK_DECOMPOSITION_H
#define CMR_BLOCK_DECOMPOSITION_H

#include <cmr/env.h>
#include "matrix_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Information on one block of a block decomposition of a matrix.
 */

typedef struct
{
  CMR_MATRIX* matrix;         /**< \brief Sparse matrix. */
  CMR_MATRIX* transpose;      /**< \brief Sparse transposed matrix. */
  size_t* rowsToOriginal;     /**< \brief Maps component rows to original matrix rows. */
  size_t* columnsToOriginal;  /**< \brief Maps component columns to original matrix columns. */
} CMR_BLOCK;

/**
 * \brief Decomposes int matrix into 1-connected submatrices.
 *
 * Allocates an array \p components with an entry per 1-connected submatrix. The caller has to free this array and
 * its members.
 */

CMR_ERROR CMRdecomposeBlocks(
  CMR* cmr,                     /**< \ref CMR environment */
  CMR_MATRIX* matrix,           /**< Matrix */
  size_t matrixType,            /**< Size of base type of matrix. */
  size_t targetType,            /**< Size of base type of component matrices. */
  size_t* pnumBlocks,           /**< Pointer for storing the number of components. */
  CMR_BLOCK** pblocks,          /**< Pointer for storing the array with component information. */
  size_t* rowsToBlock,          /**< Mapping of rows of \p matrix to components (may be \c NULL). */
  size_t* columnsToBlock,       /**< Mapping of columns of \p matrix to components (may be \c NULL). */
  size_t* rowsToBlockRows,      /**< Mapping of rows to rows of the component (may be \c NULL). */
  size_t* columnsToBlockColumns /**< Mapping of columns to columns of the component (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_BLOCK_DECOMPOSITION_H */

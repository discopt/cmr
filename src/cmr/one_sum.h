#ifndef CMR_ONESUM_INTERNAL_H
#define CMR_ONESUM_INTERNAL_H

#include <cmr/env.h>
#include "matrix_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Information on component of 1-sum of matrices.
 */

typedef struct
{
  CMR_MATRIX* matrix;      /**< \brief Sparse matrix. */
  CMR_MATRIX* transpose;   /**< \brief Sparse transposed matrix. */
  int* rowsToOriginal;    /**< \brief Maps component rows to original matrix rows. */
  int* columnsToOriginal; /**< \brief Maps component columns to original matrix columns. */
} CMR_ONESUM_COMPONENT;

/**
 * \brief Decomposes int matrix into 1-connected submatrices.
 */

CMR_ERROR decomposeOneSum(
  CMR* cmr,                           /**< \ref CMR environment */
  CMR_MATRIX* matrix,                /**< Matrix */
  size_t matrixType,                /**< Size of base type of matrix. */
  size_t targetType,                /**< Size of base type of component matrices. */
  int* pnumComponents,               /**< Number of components */
  CMR_ONESUM_COMPONENT** components, /**< Component information */
  int* rowsToComponents,            /**< Mapping of rows of \p matrix to components. Can be \c NULL. */
  int* columnsToComponents,         /**< Mapping of columns of \p matrix to components. Can be \c NULL. */
  int* rowsToComponentRows,         /**< Mapping of rows to rows of the component. Can be \c NULL. */
  int* columnsToComponentColumns    /**< Mapping of columns to columns of the component. Can be \c NULL. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_ONESUM_INTERNAL_H */

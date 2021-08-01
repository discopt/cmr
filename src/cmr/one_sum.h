#ifndef TU_ONESUM_INTERNAL_H
#define TU_ONESUM_INTERNAL_H

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
  TU_MATRIX* matrix;      /**< \brief Sparse matrix. */
  TU_MATRIX* transpose;   /**< \brief Sparse transposed matrix. */
  int* rowsToOriginal;    /**< \brief Maps component rows to original matrix rows. */
  int* columnsToOriginal; /**< \brief Maps component columns to original matrix columns. */
} TU_ONESUM_COMPONENT;

/**
 * \brief Decomposes int matrix into 1-connected submatrices.
 */

CMR_ERROR decomposeOneSum(
  TU* tu,                           /**< \ref TU environment */
  TU_MATRIX* matrix,                /**< Matrix */
  size_t matrixType,                /**< Size of base type of matrix. */
  size_t targetType,                /**< Size of base type of component matrices. */
  int* pnumComponents,               /**< Number of components */
  TU_ONESUM_COMPONENT** components, /**< Component information */
  int* rowsToComponents,            /**< Mapping of rows of \p matrix to components. Can be \c NULL. */
  int* columnsToComponents,         /**< Mapping of columns of \p matrix to components. Can be \c NULL. */
  int* rowsToComponentRows,         /**< Mapping of rows to rows of the component. Can be \c NULL. */
  int* columnsToComponentColumns    /**< Mapping of columns to columns of the component. Can be \c NULL. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_ONESUM_INTERNAL_H */

#ifndef TU_ONESUM_INTERNAL_H
#define TU_ONESUM_INTERNAL_H

#include <tu/env.h>
#include "matrix_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Information on component of 1-sum of matrices.
 */

typedef struct
{
  TU_MATRIX* matrix;          /**< \brief Sparse matrix. */
  TU_MATRIX* transpose;       /**< \brief Sparse transposed matrix. */
  size_t* rowsToOriginal;     /**< \brief Maps component rows to original matrix rows. */
  size_t* columnsToOriginal;  /**< \brief Maps component columns to original matrix columns. */
} TU_ONESUM_COMPONENT;

/**
 * \brief Decomposes int matrix into 1-connected submatrices.
 *
 * Allocates an array \p components with an entry per 1-connected submatrix. The caller has to free this array and
 * its members.
 */

TU_ERROR decomposeOneSum(
  TU* tu,                           /**< \ref TU environment */
  TU_MATRIX* matrix,                /**< Matrix */
  size_t matrixType,                /**< Size of base type of matrix. */
  size_t targetType,                /**< Size of base type of component matrices. */
  size_t* pnumComponents,           /**< Pointer for storing the number of components. */
  TU_ONESUM_COMPONENT** components, /**< Pointer for storing the array with component information. */
  size_t* rowsToComponents,         /**< Mapping of rows of \p matrix to components (may be \c NULL). */
  size_t* columnsToComponents,      /**< Mapping of columns of \p matrix to components (may be \c NULL). */
  size_t* rowsToComponentRows,      /**< Mapping of rows to rows of the component (may be \c NULL). */
  size_t* columnsToComponentColumns /**< Mapping of columns to columns of the component (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_ONESUM_INTERNAL_H */

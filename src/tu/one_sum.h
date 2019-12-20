#ifndef TU_ONESUM_INTERNAL_H
#define TU_ONESUM_INTERNAL_H

#include <tu/env.h>
#include <tu/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Information on component of 1-sum of int submatrices.
 */

typedef struct
{
  TU_SPARSE_INT matrix ;/**< Sparse matrix */
  TU_SPARSE_INT transpose; /**< Sparse transposed matrix */
  int* rowsToOriginal;
  int* columnsToOriginal;
} TU_ONESUM_COMPONENT_INT;

/**
 * \brief Decomposes int matrix into 1-connected int submatrices.
 */

void decomposeOneSumIntToInt(
  TU* tu,                               /**< TU environment */
  TU_SPARSE_INT* matrix,                /**< Sparse matrix */
  int* numComponents,                   /**< Number of components */
  TU_ONESUM_COMPONENT_INT** components, /**< Component information */
  int* rowsToComponents,                /**< Mapping of rows of \p matrix to components. Can be \c NULL. */
  int* columnsToComponents,             /**< Mapping of columns of \p matrix to components. Can be \c NULL. */
  int* rowsToComponentRows,             /**< Mapping of rows to rows of the component. Can be \c NULL. */
  int* columnsToComponentColumns        /**< Mapping of columns to columns of the component. Can be \c NULL. */
);

/**
 * \brief Information on component of 1-sum of char submatrices.
 */

typedef struct
{
  TU_SPARSE_CHAR matrix ;/**< Sparse matrix */
  TU_SPARSE_CHAR transpose; /**< Sparse transposed matrix */
  int* rowsToOriginal;
  int* columnsToOriginal;
} TU_ONESUM_COMPONENT_CHAR;

/**
 * \brief Decomposes char matrix into 1-connected char submatrices.
 */

void decomposeOneSumCharToChar(
  TU* tu,                                 /**< TU environment */
  TU_SPARSE_CHAR* matrix,                 /**< Sparse matrix */
  int* numComponents,                     /**< Number of components */
  TU_ONESUM_COMPONENT_CHAR** components,  /**< Component information */
  int* rowsToComponents,                  /**< Mapping of rows of \p matrix to components. Can be \c NULL. */
  int* columnsToComponents,               /**< Mapping of columns of \p matrix to components. Can be \c NULL. */
  int* rowsToComponentRows,               /**< Mapping of rows to rows of the component. Can be \c NULL. */
  int* columnsToComponentColumns          /**< Mapping of columns to columns of the component. Can be \c NULL. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_ONESUM_INTERNAL_H */

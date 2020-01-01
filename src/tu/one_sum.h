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
  TU_SPARSE_INT matrix;     /**< Sparse int matrix. */
  TU_SPARSE_INT transpose;  /**< Sparse transposed int matrix. */
  int* rowsToOriginal;      /**< Maps component rows to original matrix rows. */
  int* columnsToOriginal;   /**< Maps component columns to original matrix columns. */
} TU_ONESUM_COMPONENT_INT;

/**
 * \brief Decomposes int matrix into 1-connected int submatrices.
 */

void decomposeOneSumIntToInt(
  TU* tu,                               /**< TU environment */
  TU_SPARSE_INT* matrix,                /**< Sparse int matrix */
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
  TU_SPARSE_CHAR matrix;    /**< Sparse char matrix. */
  TU_SPARSE_CHAR transpose; /**< Sparse transposed char matrix. */
  int* rowsToOriginal;      /**< Maps component rows to original matrix rows. */
  int* columnsToOriginal;   /**< Maps component columns to original matrix columns. */
} TU_ONESUM_COMPONENT_CHAR;

/**
 * \brief Decomposes double matrix into 1-connected char submatrices.
 */

void decomposeOneSumDoubleToChar(
  TU* tu,                                 /**< TU environment */
  TU_SPARSE_DOUBLE* matrix,               /**< Sparse double matrix */
  int* numComponents,                     /**< Number of components */
  TU_ONESUM_COMPONENT_CHAR** components,  /**< Component information */
  int* rowsToComponents,                  /**< Mapping of rows of \p matrix to components. Can be \c NULL. */
  int* columnsToComponents,               /**< Mapping of columns of \p matrix to components. Can be \c NULL. */
  int* rowsToComponentRows,               /**< Mapping of rows to rows of the component. Can be \c NULL. */
  int* columnsToComponentColumns          /**< Mapping of columns to columns of the component. Can be \c NULL. */
);

/**
 * \brief Decomposes int matrix into 1-connected char submatrices.
 */

void decomposeOneSumIntToChar(
  TU* tu,                                 /**< TU environment */
  TU_SPARSE_INT* matrix,                  /**< Sparse int matrix */
  int* numComponents,                     /**< Number of components */
  TU_ONESUM_COMPONENT_CHAR** components,  /**< Component information */
  int* rowsToComponents,                  /**< Mapping of rows of \p matrix to components. Can be \c NULL. */
  int* columnsToComponents,               /**< Mapping of columns of \p matrix to components. Can be \c NULL. */
  int* rowsToComponentRows,               /**< Mapping of rows to rows of the component. Can be \c NULL. */
  int* columnsToComponentColumns          /**< Mapping of columns to columns of the component. Can be \c NULL. */
);


/**
 * \brief Decomposes char matrix into 1-connected char submatrices.
 */

void decomposeOneSumCharToChar(
  TU* tu,                                 /**< TU environment */
  TU_SPARSE_CHAR* matrix,                 /**< Sparse char matrix */
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

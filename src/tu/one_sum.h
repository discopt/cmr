#ifndef TU_ONESUM_INTERNAL_H
#define TU_ONESUM_INTERNAL_H

#include <tu/env.h>
#include <tu/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif
  
/**
 * \brief Decomposes int matrix into 1-connected int submatrices.
 */

void decomposeOneSumIntToInt(
  TU* tu,                         /**< TU environment */
  TU_SPARSE_INT* matrix,          /**< Sparse matrix */
  int* numComponents,             /**< Number of components */
  TU_SPARSE_INT** compMatrices,   /**< Array of sparse matrices of components */
  TU_SPARSE_INT** compTransposes, /**< Array of sparse transposed matrices of components */
  int*** rowMapping,              /**< Array mapping component rows to original rows */
  int*** columnMapping            /**< Array mapping component columns to original columns */
);

/**
 * \brief Decomposes char matrix into 1-connected char submatrices.
 */

void decomposeOneSumCharToChar(
  TU* tu,                           /**< TU environment */
  TU_SPARSE_CHAR* matrix,           /**< Sparse matrix */
  int* numComponents,               /**< Number of components */
  TU_SPARSE_CHAR** compMatrices,    /**< Array of sparse matrices of components */
  TU_SPARSE_CHAR** compTransposes,  /**< Array of sparse transposed matrices of components */
  int*** rowMapping,                /**< Array mapping component rows to original rows */
  int*** columnMapping              /**< Array mapping component columns to original columns */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_ONESUM_INTERNAL_H */

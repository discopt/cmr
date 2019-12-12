#ifndef TU_ONESUM_INTERNAL_H
#define TU_ONESUM_INTERNAL_H

#include <tu/env.h>
#include <tu/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

void decomposeOneSumIntInt(
  TU* tu,                         /**< TU environment */
  TU_SPARSE_INT* matrix,          /**< Sparse matrix */
  int* numComponents,             /**< Number of components */
  TU_SPARSE_INT** compMatrices,   /**< Array of sparse matrices of components */
  TU_SPARSE_INT** compTransposes, /**< Array of sparse transposed matrices of components */
  int*** rowMapping,              /**< Array mapping component rows to original rows */
  int*** columnMapping            /**< Array mapping component columns to original columns */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_ONESUM_INTERNAL_H */

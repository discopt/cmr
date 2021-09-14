#ifndef CMR_DENSEMATRIX_INTERNAL_H
#define CMR_DENSEMATRIX_INTERNAL_H

#include <cmr/env.h>
#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Dense matrix.
 */

typedef struct
{
  unsigned long long* data;
  size_t numRows;
  size_t numColumns;
} DenseBinaryMatrix;

CMR_ERROR CMRdensebinmatrixCreateStack(
  CMR* cmr,                   /**< \ref CMR environment. */
  size_t numRows,             /**< Number of rows. */
  size_t numColumns,          /**< Number of columns. */
  DenseBinaryMatrix** presult /**< Pointer for storing the result. */
);

CMR_ERROR CMRdensebinmatrixFreeStack(
  CMR* cmr,                   /**< \ref CMR environment. */
  DenseBinaryMatrix** pmatrix /**< Pointer for storing the result. */
);

static inline
bool CMRdensebinmatrixGet(
  DenseBinaryMatrix* matrix,  /**< Matrix. */
  size_t row,                 /**< Row index. */
  size_t column               /**< Column index. */
)
{
  size_t index = row * matrix->numColumns + column;
  unsigned long long block = matrix->data[index / (8 * sizeof(unsigned long long))];
  return block & (1UL << (index % (8 * sizeof(unsigned long long))));
}

static inline
void CMRdensebinmatrixSet0(
  DenseBinaryMatrix* matrix,  /**< Matrix. */
  size_t row,                 /**< Row index. */
  size_t column               /**< Column index. */
)
{
  size_t index = row * matrix->numColumns + column;
  unsigned long long* pblock = &matrix->data[index / (8 * sizeof(unsigned long long))];
  *pblock &= ~(1UL << (index % (8 * sizeof(unsigned long long))));
}

static inline
void CMRdensebinmatrixSet1(
  DenseBinaryMatrix* matrix,  /**< Matrix. */
  size_t row,                 /**< Row index. */
  size_t column               /**< Column index. */
)
{
  size_t index = row * matrix->numColumns + column;
  unsigned long long* pblock = &matrix->data[index / (8 * sizeof(unsigned long long))];
  *pblock |= (1UL << (index % (8 * sizeof(unsigned long long))));
}

static inline
void CMRdensebinmatrixSet(
  DenseBinaryMatrix* matrix,  /**< Matrix. */
  size_t row,                 /**< Row index. */
  size_t column,              /**< Column index. */
  bool value                  /**< Value. */
)
{
  size_t index = row * matrix->numColumns + column;
  unsigned long long* pblock = &matrix->data[index / (8 * sizeof(unsigned long long))];
  size_t mask = (1UL << (index % (8 * sizeof(unsigned long long))));
  if (value)
    *pblock |= mask;
  else
    *pblock &= ~mask;
}

#ifdef __cplusplus
}
#endif

#endif /* CMR_LISTMATRIX_INTERNAL_H */

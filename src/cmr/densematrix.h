#define CMR_DEBUG

#ifndef CMR_DENSEMATRIX_INTERNAL_H
#define CMR_DENSEMATRIX_INTERNAL_H

#include "env_internal.h"

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
  CMRdbgMsg(8, "CMRdensebinmatrixGet(%zu,%zu) uses index %zu in block %zu at %zu\n", row, column, index,
    index / (8 * sizeof(unsigned long long)), (index % (8 * sizeof(unsigned long long))));
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
  CMRdbgMsg(8, "CMRdensebinmatrixSet0(%zu,%zu) uses index %zu in block %zu at %zu.\n", row, column, index,
    index / (8 * sizeof(unsigned long long)), (index % (8 * sizeof(unsigned long long))));
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
  CMRdbgMsg(8, "CMRdensebinmatrixSet1(%zu,%zu) uses index %zu in block %zu at %zu.\n", row, column, index,
    index / (8 * sizeof(unsigned long long)), (index % (8 * sizeof(unsigned long long))));
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
  CMRdbgMsg(8, "CMRdensebinmatrixSet(%zu,%zu,%d) uses index %zu in block %zu at %zu.\n", row, column, value ? 1 : 0, index,
    index / (8 * sizeof(unsigned long long)), (index % (8 * sizeof(unsigned long long))));
  unsigned long long* pblock = &matrix->data[index / (8 * sizeof(unsigned long long))];
  size_t mask = (1UL << (index % (8 * sizeof(unsigned long long))));
  if (value)
    *pblock |= mask;
  else
    *pblock &= ~mask;
}


static inline
void CMRdensebinmatrixFlip(
  DenseBinaryMatrix* matrix,  /**< Matrix. */
  size_t row,                 /**< Row index. */
  size_t column               /**< Column index. */
)
{
  size_t index = row * matrix->numColumns + column;
  CMRdbgMsg(8, "CMRdensebinmatrixFlip(%zu,%zu) uses index %zu in block %zu at %zu.\n", row, column, index,
    index / (8 * sizeof(unsigned long long)), (index % (8 * sizeof(unsigned long long))));
  unsigned long long* pblock = &matrix->data[index / (8 * sizeof(unsigned long long))];
  size_t mask = (1UL << (index % (8 * sizeof(unsigned long long))));
  if (*pblock & mask)
    *pblock &= ~mask;
  else
    *pblock |= mask;
}

#ifdef __cplusplus
}
#endif

#endif /* CMR_LISTMATRIX_INTERNAL_H */

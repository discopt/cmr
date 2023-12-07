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
  size_t* data;
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
  size_t block = matrix->data[index / (8 * sizeof(size_t))];
  CMRdbgMsg(8, "CMRdensebinmatrixGet(%zu,%zu) uses index %zu in block %zu at %zu -> %d\n", row, column, index,
    index / (8 * sizeof(size_t)), (index % (8 * sizeof(size_t))),
    (block & (1UL << (index % (8 * sizeof(size_t))))) ? 1 : 0);
  return block & (1UL << (index % (8 * sizeof(size_t))));
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
    index / (8 * sizeof(size_t)), (index % (8 * sizeof(size_t))));
  size_t* pblock = &matrix->data[index / (8 * sizeof(size_t))];
  *pblock &= ~(1UL << (index % (8 * sizeof(size_t))));
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
    index / (8 * sizeof(size_t)), (index % (8 * sizeof(size_t))));
  size_t* pblock = &matrix->data[index / (8 * sizeof(size_t))];
  *pblock |= (1UL << (index % (8 * sizeof(size_t))));
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
    index / (8 * sizeof(size_t)), (index % (8 * sizeof(size_t))));
  size_t* pblock = &matrix->data[index / (8 * sizeof(size_t))];
  size_t mask = (1UL << (index % (8 * sizeof(size_t))));
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
    index / (8 * sizeof(size_t)), (index % (8 * sizeof(size_t))));
  size_t* pblock = &matrix->data[index / (8 * sizeof(size_t))];
  size_t mask = (1UL << (index % (8 * sizeof(size_t))));
  if (*pblock & mask)
    *pblock &= ~mask;
  else
    *pblock |= mask;
}

#ifdef __cplusplus
}
#endif

#endif /* CMR_LISTMATRIX_INTERNAL_H */

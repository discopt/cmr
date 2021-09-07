#ifndef CMR_LIST_MATRIX_INTERNAL_H
#define CMR_LIST_MATRIX_INTERNAL_H

#include <cmr/env.h>
#include <cmr/matrix.h>

#include "hashtable.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Nonzero of \ref ListMatrix.
 */

typedef struct _ListMatrixNonzero
{
  struct _ListMatrixNonzero* left;  /**< \brief Pointer to previous nonzero in the same row. */
  struct _ListMatrixNonzero* right; /**< \brief Pointer to next nonzero in the same row. */
  struct _ListMatrixNonzero* above; /**< \brief Pointer to previous nonzero in the same column. */
  struct _ListMatrixNonzero* below; /**< \brief Pointer to next nonzero in the same column. */
  size_t row;                       /**< \brief Row. */
  size_t column;                    /**< \brief Column. */
  char value;                       /**< \brief Matrix entry. */
  long special : 56;                /**< \brief Remaining bits (on 64 bit) may be used for a special purpose. */
} ListMatrixNonzero;

/**
 * \brief Row/column information of \ref ListMatrix.
 */

typedef struct
{
  ListMatrixNonzero head;             /**< \brief Dummy nonzero in that row/column. */
  size_t numNonzeros;                 /**< \brief Number of nonzeros in that row/column. */
} ListMatrixElement;

/**
 * \brief Linked-list representation of a matrix.
 *
 * The Each nonzero is part of two doubly-linked lists, one for all nonzeros in the same row and one for all the
 * nonzeros in the same column.
 *
 * \note If the allocated memory for the rows, columns or nonzeros does not suffice, it is automatically reallocated.
 *       However, this means that a pointer to a \ref ListMatrixNonzero struct is only valid if no reallocation occurs.
 *       Hence, the user must make sure that sufficient memory is allocated or that no such pointers are used.
 */

typedef struct
{
  size_t memRows;                       /**< \brief Memory for rows. */
  size_t numRows;                       /**< \brief Number of rows. */
  ListMatrixElement* rowElements;       /**< \brief Row data. */
  size_t memColumns;                    /**< \brief Memory for columns. */
  size_t numColumns;                    /**< \brief Number of columns. */
  ListMatrixElement* columnElements;    /**< \brief Column data. */

  size_t numNonzeros;
  ListMatrixNonzero anchor;             /**< \brief Anchor for nonzeros. */
  size_t memNonzeros;                   /**< \brief Amount of memory for nonzeros. */
  ListMatrixNonzero* nonzeros;          /**< \brief Raw nonzero data. */
  ListMatrixNonzero* firstFreeNonzero;  /**< \brief Beginning of free list. */
} ListMatrix;

/**
 * \brief Allocates memory for a list matrix.
 */

CMR_ERROR CMRlistmatrixAlloc(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t memRows,       /**< Memory for rows. */
  size_t memColumns,    /**< Memory for columns. */
  size_t memNonzeros,   /**< Memory for nonzeros. */
  ListMatrix** presult  /**< Pointer for storing the created list matrix. */
);

/**
 * \brief Frees a list matrix.
 */

CMR_ERROR CMRlistmatrixFree(
  CMR* cmr,                 /**< \ref CMR environment. */
  ListMatrix** plistmatrix  /**< Pointer to list matrix. */
);

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmatrixInitializeFromMatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatrix* listmatrix, /**< List matrix. */
  CMR_CHRMAT* matrix      /**< Matrix to be copied to \p listmatrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_LIST_MATRIX_INTERNAL_H */

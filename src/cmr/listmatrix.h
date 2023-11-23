#ifndef CMR_LISTMATRIX_INTERNAL_H
#define CMR_LISTMATRIX_INTERNAL_H

#include <cmr/env.h>
#include <cmr/matrix.h>

#include "hashtable.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Nonzero of a \ref ChrListMat.
 */

typedef struct _ChrListMatNonzero
{
  struct _ChrListMatNonzero* left;  /**< \brief Pointer to previous nonzero in the same row. */
  struct _ChrListMatNonzero* right; /**< \brief Pointer to next nonzero in the same row. */
  struct _ChrListMatNonzero* above; /**< \brief Pointer to previous nonzero in the same column. */
  struct _ChrListMatNonzero* below; /**< \brief Pointer to next nonzero in the same column. */
  size_t row;                       /**< \brief Row. */
  size_t column;                    /**< \brief Column. */
  char value;                       /**< \brief Matrix entry. */
  long special : 56;                /**< \brief Remaining bits (on 64 bit) may be used for a special purpose. */
} ChrListMatNonzero;

/**
 * \brief Nonzero of a \ref LongListMat.
 */

typedef struct _LongListMatNonzero
{
  struct _LongListMatNonzero* left;   /**< \brief Pointer to previous nonzero in the same row. */
  struct _LongListMatNonzero* right;  /**< \brief Pointer to next nonzero in the same row. */
  struct _LongListMatNonzero* above;  /**< \brief Pointer to previous nonzero in the same column. */
  struct _LongListMatNonzero* below;  /**< \brief Pointer to next nonzero in the same column. */
  size_t row;                         /**< \brief Row. */
  size_t column;                      /**< \brief Column. */
  int64_t value;                      /**< \brief Matrix entry. */
  long special : 56;                  /**< \brief Remaining bits (on 64 bit) may be used for a special purpose. */
} LongListMatNonzero;

/**
 * \brief Row/column information of a \ref ChrListMat.
 */

typedef struct
{
  ChrListMatNonzero head; /**< \brief Dummy nonzero in that row/column. */
  size_t numNonzeros;     /**< \brief Number of nonzeros in that row/column. */
} ChrListMatElement;

/**
 * \brief Row/column information of a \ref LongListMat.
 */

typedef struct
{
  LongListMatNonzero head;  /**< \brief Dummy nonzero in that row/column. */
  size_t numNonzeros;       /**< \brief Number of nonzeros in that row/column. */
} LongListMatElement;

/**
 * \brief Linked-list representation of a matrix with char values.
 *
 * The Each nonzero is part of two doubly-linked lists, one for all nonzeros in the same row and one for all the
 * nonzeros in the same column.
 *
 * \note If the allocated memory for the rows, columns or nonzeros does not suffice, it is automatically reallocated.
 *       However, this means that a pointer to a \ref ChrListMatNonzero struct is only valid if no reallocation occurs.
 *       Hence, the user must make sure that sufficient memory is allocated or that no such pointers are used.
 */

typedef struct
{
  size_t memRows;                       /**< \brief Memory for rows. */
  size_t numRows;                       /**< \brief Number of rows. */
  ChrListMatElement* rowElements;       /**< \brief Row data. */
  size_t memColumns;                    /**< \brief Memory for columns. */
  size_t numColumns;                    /**< \brief Number of columns. */
  ChrListMatElement* columnElements;    /**< \brief Column data. */

  size_t numNonzeros;
  ChrListMatNonzero anchor;             /**< \brief Anchor for nonzeros. */
  size_t memNonzeros;                   /**< \brief Amount of memory for nonzeros. */
  ChrListMatNonzero* nonzeros;          /**< \brief Raw nonzero data. */
  ChrListMatNonzero* firstFreeNonzero;  /**< \brief Beginning of free list. */
} ChrListMat;

/**
 * \brief Linked-list representation of a matrix with 64-bit integer values.
 *
 * The Each nonzero is part of two doubly-linked lists, one for all nonzeros in the same row and one for all the
 * nonzeros in the same column.
 *
 * \note If the allocated memory for the rows, columns or nonzeros does not suffice, it is automatically reallocated.
 *       However, this means that a pointer to a \ref LongListMatNonzero struct is only valid if no reallocation occurs.
 *       Hence, the user must make sure that sufficient memory is allocated or that no such pointers are used.
 */

typedef struct
{
  size_t memRows;                       /**< \brief Memory for rows. */
  size_t numRows;                       /**< \brief Number of rows. */
  LongListMatElement* rowElements;      /**< \brief Row data. */
  size_t memColumns;                    /**< \brief Memory for columns. */
  size_t numColumns;                    /**< \brief Number of columns. */
  LongListMatElement* columnElements;   /**< \brief Column data. */

  size_t numNonzeros;
  LongListMatNonzero anchor;            /**< \brief Anchor for nonzeros. */
  size_t memNonzeros;                   /**< \brief Amount of memory for nonzeros. */
  LongListMatNonzero* nonzeros;         /**< \brief Raw nonzero data. */
  LongListMatNonzero* firstFreeNonzero; /**< \brief Beginning of free list. */
} LongListMat;

/**
 * \brief Allocates memory for a chr list matrix.
 */

CMR_ERROR CMRchrlistmatAlloc(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t memRows,       /**< Memory for rows. */
  size_t memColumns,    /**< Memory for columns. */
  size_t memNonzeros,   /**< Memory for nonzeros. */
  ChrListMat** presult  /**< Pointer for storing the created list matrix. */
);

/**
 * \brief Allocates memory for a long list matrix.
 */

CMR_ERROR CMRlonglistmatAlloc(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t memRows,       /**< Memory for rows. */
  size_t memColumns,    /**< Memory for columns. */
  size_t memNonzeros,   /**< Memory for nonzeros. */
  LongListMat** presult /**< Pointer for storing the created list matrix. */
);

/**
 * \brief Frees a chr list matrix.
 */

CMR_ERROR CMRchrlistmatFree(
  CMR* cmr,                 /**< \ref CMR environment. */
  ChrListMat** plistmatrix  /**< Pointer to list matrix. */
);

/**
 * \brief Frees a long list matrix.
 */

CMR_ERROR CMRlonglistmatFree(
  CMR* cmr,                 /**< \ref CMR environment. */
  LongListMat** plistmatrix /**< Pointer to list matrix. */
);

/**
 * \brief Initializes a zero chr list matrix.
 */

CMR_ERROR CMRchrlistmatInitializeZero(
  CMR* cmr,               /**< \ref CMR environment. */
  ChrListMat* listmatrix, /**< List matrix. */
  size_t numRows,         /**< Number of rows. */
  size_t numColumns       /**< Number of columns. */
);

/**
 * \brief Initializes a zero long list matrix.
 */

CMR_ERROR CMRlonglistmatInitializeZero(
  CMR* cmr,                 /**< \ref CMR environment. */
  LongListMat* listmatrix,  /**< List matrix. */
  size_t numRows,           /**< Number of rows. */
  size_t numColumns         /**< Number of columns. */
);

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRchrlistmatInitializeFromChrMatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ChrListMat* listmatrix, /**< List matrix. */
  CMR_CHRMAT* matrix      /**< Matrix to be copied to \p listmatrix. */
);

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlonglistmatInitializeFromIntMatrix(
  CMR* cmr,                 /**< \ref CMR environment. */
  LongListMat* listmatrix,  /**< List matrix. */
  CMR_INTMAT* matrix        /**< Matrix to be copied to \p listmatrix. */
);

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRchrlistmatInitializeFromDoubleMatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ChrListMat* listmatrix, /**< List matrix. */
  CMR_DBLMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  double epsilon          /**< Tolerance to consider as exact integer. */
);

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlonglistmatInitializeFromDoubleMatrix(
  CMR* cmr,                 /**< \ref CMR environment. */
  LongListMat* listmatrix,  /**< List matrix. */
  CMR_DBLMAT* matrix,       /**< Matrix to be copied to \p listmatrix. */
  double epsilon            /**< Tolerance to consider as exact integer. */
);

/**
 * \brief Copies \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRchrlistmatInitializeFromChrSubmatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ChrListMat* listmatrix, /**< List matrix. */
  CMR_CHRMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix   /**< Submatrix of \p matrix. */
);

/**
 * \brief Copies \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlonglistmatInitializeFromIntSubmatrix(
  CMR* cmr,                 /**< \ref CMR environment. */
  LongListMat* listmatrix,  /**< List matrix. */
  CMR_INTMAT* matrix,       /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix     /**< Submatrix of \p matrix. */
);

/**
 * \brief Copies all but \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRchrlistmatInitializeFromSubmatrixComplement(
  CMR* cmr,               /**< \ref CMR environment. */
  ChrListMat* listmatrix, /**< List matrix. */
  CMR_CHRMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix   /**< Submatrix of \p matrix. */
);

/**
 * \brief Copies all but \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlonglistmatInitializeFromIntSubmatrixComplement(
  CMR* cmr,                 /**< \ref CMR environment. */
  LongListMat* listmatrix,  /**< List matrix. */
  CMR_INTMAT* matrix,       /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix     /**< Submatrix of \p matrix. */
);

/**
 * \brief Prints the chr list matrix as a dense matrix.
 */

CMR_ERROR CMRchrlistmatPrintDense(
  CMR* cmr,               /**< \ref CMR environment. */
  ChrListMat* listmatrix, /**< List matrix. */
  FILE* stream            /**< Stream to print to. */
);

/**
 * \brief Prints the long list matrix as a dense matrix.
 */

CMR_ERROR CMRlonglistmatPrintDense(
  CMR* cmr,                 /**< \ref CMR environment. */
  LongListMat* listmatrix,  /**< List matrix. */
  FILE* stream              /**< Stream to print to. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_LISTMATRIX_INTERNAL_H */

#ifndef CMR_LISTMATRIX_INTERNAL_H
#define CMR_LISTMATRIX_INTERNAL_H

#include <stddef.h>
#include <inttypes.h>

#include <cmr/env.h>
#include <cmr/matrix.h>

#include "hashtable.h"

#if defined(CMR_WITH_GMP)
#include <gmp.h>
#endif /* CMR_WITH_GMP */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Nonzero of a \ref ListMat8.
 */

typedef struct _ListMat8Nonzero
{
  struct _ListMat8Nonzero* left;  /**< \brief Pointer to previous nonzero in the same row. */
  struct _ListMat8Nonzero* right; /**< \brief Pointer to next nonzero in the same row. */
  struct _ListMat8Nonzero* above; /**< \brief Pointer to previous nonzero in the same column. */
  struct _ListMat8Nonzero* below; /**< \brief Pointer to next nonzero in the same column. */
  size_t row;                     /**< \brief Row. */
  size_t column;                  /**< \brief Column. */
  int8_t value;                   /**< \brief Matrix entry. */
  int64_t special : 56;           /**< \brief Remaining bits (on 64 bit) may be used for a special purpose. */
} ListMat8Nonzero;

/**
 * \brief Row/column information of a \ref ListMat8.
 */

typedef struct
{
  ListMat8Nonzero head; /**< \brief Dummy nonzero in that row/column. */
  size_t numNonzeros;   /**< \brief Number of nonzeros in that row/column. */
} ListMat8Element;

/**
 * \brief Linked-list representation of a matrix with 8-bit integer values.
 *
 * The Each nonzero is part of two doubly-linked lists, one for all nonzeros in the same row and one for all the
 * nonzeros in the same column.
 *
 * \note If the allocated memory for the rows, columns or nonzeros does not suffice, it is automatically reallocated.
 *       However, this means that a pointer to a \ref ListMat8Nonzero struct is only valid if no reallocation occurs.
 *       Hence, the user must make sure that sufficient memory is allocated or that no such pointers are used.
 */

typedef struct
{
  size_t memRows;                     /**< \brief Memory for rows. */
  size_t numRows;                     /**< \brief Number of rows. */
  ListMat8Element* rowElements;       /**< \brief Row data. */
  size_t memColumns;                  /**< \brief Memory for columns. */
  size_t numColumns;                  /**< \brief Number of columns. */
  ListMat8Element* columnElements;    /**< \brief Column data. */

  size_t numNonzeros;
  ListMat8Nonzero anchor;             /**< \brief Anchor for nonzeros. */
  size_t memNonzeros;                 /**< \brief Amount of memory for nonzeros. */
  ListMat8Nonzero* nonzeros;          /**< \brief Raw nonzero data. */
  ListMat8Nonzero* firstFreeNonzero;  /**< \brief Beginning of free list; uses \c right pointers. */
} ListMat8;


/**
 * \brief Nonzero of a \ref ListMat64.
 */

typedef struct _ListMat64Nonzero
{
  struct _ListMat64Nonzero* left;   /**< \brief Pointer to previous nonzero in the same row. */
  struct _ListMat64Nonzero* right;  /**< \brief Pointer to next nonzero in the same row. */
  struct _ListMat64Nonzero* above;  /**< \brief Pointer to previous nonzero in the same column. */
  struct _ListMat64Nonzero* below;  /**< \brief Pointer to next nonzero in the same column. */
  size_t row;                       /**< \brief Row. */
  size_t column;                    /**< \brief Column. */
  int64_t value;                    /**< \brief Matrix entry. */
  long special;                     /**< \brief May be used for a special purpose. */
} ListMat64Nonzero;

/**
 * \brief Row/column information of a \ref ListMat64.
 */

typedef struct
{
  ListMat64Nonzero head;  /**< \brief Dummy nonzero in that row/column. */
  size_t numNonzeros;     /**< \brief Number of nonzeros in that row/column. */
} ListMat64Element;

/**
 * \brief Linked-list representation of a matrix with 64-bit integer values.
 *
 * The Each nonzero is part of two doubly-linked lists, one for all nonzeros in the same row and one for all the
 * nonzeros in the same column.
 *
 * \note If the allocated memory for the rows, columns or nonzeros does not suffice, it is automatically reallocated.
 *       However, this means that a pointer to a \ref ListMat64Nonzero struct is only valid if no reallocation occurs.
 *       Hence, the user must make sure that sufficient memory is allocated or that no such pointers are used.
 */

typedef struct
{
  size_t memRows;                     /**< \brief Memory for rows. */
  size_t numRows;                     /**< \brief Number of rows. */
  ListMat64Element* rowElements;      /**< \brief Row data. */
  size_t memColumns;                  /**< \brief Memory for columns. */
  size_t numColumns;                  /**< \brief Number of columns. */
  ListMat64Element* columnElements;   /**< \brief Column data. */

  size_t numNonzeros;
  ListMat64Nonzero anchor;            /**< \brief Anchor for nonzeros. */
  size_t memNonzeros;                 /**< \brief Amount of memory for nonzeros. */
  ListMat64Nonzero* nonzeros;         /**< \brief Raw nonzero data. */
  ListMat64Nonzero* firstFreeNonzero; /**< \brief Beginning of free list; uses \c right pointers. */
} ListMat64;

#if defined(CMR_WITH_GMP)

/**
 * \brief Nonzero of a \ref ListMat64.
 */

typedef struct _ListMatGMPNonzero
{
  struct _ListMatGMPNonzero* left;   /**< \brief Pointer to previous nonzero in the same row. */
  struct _ListMatGMPNonzero* right;  /**< \brief Pointer to next nonzero in the same row. */
  struct _ListMatGMPNonzero* above;  /**< \brief Pointer to previous nonzero in the same column. */
  struct _ListMatGMPNonzero* below;  /**< \brief Pointer to next nonzero in the same column. */
  size_t row;                       /**< \brief Row. */
  size_t column;                    /**< \brief Column. */
  mpz_t value;                      /**< \brief Matrix entry. */
  long special;                     /**< \brief May be used for a special purpose. */
} ListMatGMPNonzero;

/**
 * \brief Row/column information of a \ref ListMat64.
 */

typedef struct
{
  ListMatGMPNonzero head;  /**< \brief Dummy nonzero in that row/column. */
  size_t numNonzeros;     /**< \brief Number of nonzeros in that row/column. */
} ListMatGMPElement;

/**
 * \brief Linked-list representation of a matrix with 64-bit integer values.
 *
 * The Each nonzero is part of two doubly-linked lists, one for all nonzeros in the same row and one for all the
 * nonzeros in the same column.
 *
 * \note If the allocated memory for the rows, columns or nonzeros does not suffice, it is automatically reallocated.
 *       However, this means that a pointer to a \ref ListMat64Nonzero struct is only valid if no reallocation occurs.
 *       Hence, the user must make sure that sufficient memory is allocated or that no such pointers are used.
 */

typedef struct
{
  size_t memRows;                     /**< \brief Memory for rows. */
  size_t numRows;                     /**< \brief Number of rows. */
  ListMatGMPElement* rowElements;      /**< \brief Row data. */
  size_t memColumns;                  /**< \brief Memory for columns. */
  size_t numColumns;                  /**< \brief Number of columns. */
  ListMatGMPElement* columnElements;   /**< \brief Column data. */

  size_t numNonzeros;
  ListMatGMPNonzero anchor;            /**< \brief Anchor for nonzeros. */
  size_t memNonzeros;                 /**< \brief Amount of memory for nonzeros. */
  ListMatGMPNonzero* nonzeros;         /**< \brief Raw nonzero data. */
  ListMatGMPNonzero* firstFreeNonzero; /**< \brief Beginning of free list; uses \c right pointers. */
} ListMatGMP;

#endif /* CMR_WITH_GMP */

/**
 * \brief Allocates memory for an 8-bit list matrix.
 */

CMR_ERROR CMRlistmat8Alloc(
  CMR* cmr,           /**< \ref CMR environment. */
  size_t memRows,     /**< Memory for rows. */
  size_t memColumns,  /**< Memory for columns. */
  size_t memNonzeros, /**< Memory for nonzeros. */
  ListMat8** presult  /**< Pointer for storing the created list matrix. */
);

/**
 * \brief Allocates memory for a 64-bit list matrix.
 */

CMR_ERROR CMRlistmat64Alloc(
  CMR* cmr,           /**< \ref CMR environment. */
  size_t memRows,     /**< Memory for rows. */
  size_t memColumns,  /**< Memory for columns. */
  size_t memNonzeros, /**< Memory for nonzeros. */
  ListMat64** presult /**< Pointer for storing the created list matrix. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Allocates memory for a GMP list matrix.
 */

CMR_ERROR CMRlistmatGMPAlloc(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t memRows,       /**< Memory for rows. */
  size_t memColumns,    /**< Memory for columns. */
  size_t memNonzeros,   /**< Memory for nonzeros. */
  ListMatGMP** presult  /**< Pointer for storing the created list matrix. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Frees an 8-bit list matrix.
 */

CMR_ERROR CMRlistmat8Free(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat8** plistmatrix  /**< Pointer to list matrix. */
);

/**
 * \brief Frees a 64-bit list matrix.
 */

CMR_ERROR CMRlistmat64Free(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64** plistmatrix /**< Pointer to list matrix. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Frees a gmp list matrix.
 */

CMR_ERROR CMRlistmatGMPFree(
  CMR* cmr,                 /**< \ref CMR environment. */
  ListMatGMP** plistmatrix  /**< Pointer to list matrix. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Initializes a zero 8-bit list matrix.
 */

CMR_ERROR CMRlistmat8InitializeZero(
  CMR* cmr,             /**< \ref CMR environment. */
  ListMat8* listmatrix, /**< List matrix. */
  size_t numRows,       /**< Number of rows. */
  size_t numColumns     /**< Number of columns. */
);

/**
 * \brief Initializes a zero 64-bit list matrix.
 */

CMR_ERROR CMRlistmat64InitializeZero(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64* listmatrix,  /**< List matrix. */
  size_t numRows,         /**< Number of rows. */
  size_t numColumns       /**< Number of columns. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Initializes a zero 64-bit list matrix.
 */

CMR_ERROR CMRlistmatGMPInitializeZero(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatGMP* listmatrix, /**< List matrix. */
  size_t numRows,         /**< Number of rows. */
  size_t numColumns       /**< Number of columns. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmat8InitializeFromChrMatrix(
  CMR* cmr,             /**< \ref CMR environment. */
  ListMat8* listmatrix, /**< List matrix. */
  CMR_CHRMAT* matrix    /**< Matrix to be copied to \p listmatrix. */
);

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmat64InitializeFromIntMatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64* listmatrix,  /**< List matrix. */
  CMR_INTMAT* matrix      /**< Matrix to be copied to \p listmatrix. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmatGMPInitializeFromIntMatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatGMP* listmatrix, /**< List matrix. */
  CMR_INTMAT* matrix      /**< Matrix to be copied to \p listmatrix. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmat8InitializeFromDoubleMatrix(
  CMR* cmr,             /**< \ref CMR environment. */
  ListMat8* listmatrix, /**< List matrix. */
  CMR_DBLMAT* matrix,   /**< Matrix to be copied to \p listmatrix. */
  double epsilon        /**< Tolerance to consider as exact integer. */
);

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmat64InitializeFromDoubleMatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64* listmatrix,  /**< List matrix. */
  CMR_DBLMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  double epsilon          /**< Tolerance to consider as exact integer. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Copies \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmatGMPInitializeFromDoubleMatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatGMP* listmatrix,  /**< List matrix. */
  CMR_DBLMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  double epsilon          /**< Tolerance to consider as exact integer. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Copies \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmat8InitializeFromChrSubmatrix(
  CMR* cmr,             /**< \ref CMR environment. */
  ListMat8* listmatrix, /**< List matrix. */
  CMR_CHRMAT* matrix,   /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix /**< Submatrix of \p matrix. */
);

/**
 * \brief Copies \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmat64InitializeFromIntSubmatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64* listmatrix,  /**< List matrix. */
  CMR_INTMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix   /**< Submatrix of \p matrix. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Copies \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmatGMPInitializeFromIntSubmatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatGMP* listmatrix, /**< List matrix. */
  CMR_INTMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix   /**< Submatrix of \p matrix. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Copies all but \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmat8InitializeFromSubmatrixComplement(
  CMR* cmr,             /**< \ref CMR environment. */
  ListMat8* listmatrix, /**< List matrix. */
  CMR_CHRMAT* matrix,   /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix /**< Submatrix of \p matrix. */
);

/**
 * \brief Copies all but \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmat64InitializeFromIntSubmatrixComplement(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64* listmatrix,  /**< List matrix. */
  CMR_INTMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix   /**< Submatrix of \p matrix. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Copies all but \p submatrix of \p matrix into \p listmatrix.
 */

CMR_ERROR CMRlistmatGMPInitializeFromIntSubmatrixComplement(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatGMP* listmatrix, /**< List matrix. */
  CMR_INTMAT* matrix,     /**< Matrix to be copied to \p listmatrix. */
  CMR_SUBMAT* submatrix   /**< Submatrix of \p matrix. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Prints the 8-bit list matrix as a dense matrix.
 */

CMR_ERROR CMRlistmat8PrintDense(
  CMR* cmr,             /**< \ref CMR environment. */
  ListMat8* listmatrix, /**< List matrix. */
  FILE* stream          /**< Stream to print to. */
);

/**
 * \brief Prints the 64-bit list matrix as a dense matrix.
 */

CMR_ERROR CMRlistmat64PrintDense(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64* listmatrix,  /**< List matrix. */
  FILE* stream            /**< Stream to print to. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Prints the GMP list matrix as a dense matrix.
 */

CMR_ERROR CMRlistmatGMPPrintDense(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatGMP* listmatrix, /**< List matrix. */
  FILE* stream            /**< Stream to print to. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Creates a new element and inserts it into the doubly-linked lists.
 *
 * The function may reallocate the array of nonzeros.
 **/

CMR_ERROR CMRlistmat8Insert(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat8* listmatrix,   /**< List matrix. */
  size_t row,             /**< Row of new element. */
  size_t column,          /**< Column of new element. */
  int8_t value,           /**< Value of new element. */
  long special,           /**< Special entry of new element. */
  ptrdiff_t* pmemoryShift /**< If not \c NULL, each nonzero's address is shifted by this value. */
);

/**
 * \brief Creates a new element and inserts it into the doubly-linked lists.
 *
 * The function may reallocate the array of nonzeros.
 **/

CMR_ERROR CMRlistmat64Insert(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64* listmatrix,  /**< List matrix. */
  size_t row,             /**< Row of new element. */
  size_t column,          /**< Column of new element. */
  int64_t value,          /**< Value of new element. */
  long special,           /**< Special entry of new element. */
  ptrdiff_t* pmemoryShift /**< If not \c NULL, each nonzero's address is shifted by this value. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Creates a new element and inserts it into the doubly-linked lists.
 *
 * The function may reallocate the array of nonzeros.
 **/

CMR_ERROR CMRlistmatGMPInsert(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatGMP* listmatrix, /**< List matrix. */
  size_t row,             /**< Row of new element. */
  size_t column,          /**< Column of new element. */
  mpz_srcptr value,       /**< Value of new element. */
  long special,           /**< Special entry of new element. */
  ptrdiff_t* pmemoryShift /**< If not \c NULL, each nonzero's address is shifted by this value. */
);

#endif /* CMR_WITH_GMP */

/**
 * \brief Delete a nonzero element.
 **/

CMR_ERROR CMRlistmat8Delete(
  CMR* cmr,             /**< \ref CMR environment. */
  ListMat8* listmatrix, /**< List matrix. */
    ListMat8Nonzero* nz /**< Nonzero to delete. */
);

/**
 * \brief Delete a nonzero element.
 **/

CMR_ERROR CMRlistmat64Delete(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMat64* listmatrix,  /**< List matrix. */
  ListMat64Nonzero* nz    /**< Nonzero to delete. */
);

#if defined(CMR_WITH_GMP)

/**
 * \brief Delete a nonzero element.
 **/

CMR_ERROR CMRlistmatGMPDelete(
  CMR* cmr,               /**< \ref CMR environment. */
  ListMatGMP* listmatrix, /**< List matrix. */
  ListMatGMPNonzero* nz   /**< Nonzero to delete. */
);

#endif /* CMR_WITH_GMP */

#ifdef __cplusplus
}
#endif

#endif /* CMR_LISTMATRIX_INTERNAL_H */

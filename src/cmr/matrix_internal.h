#ifndef CMR_MATRIX_INTERNAL_H
#define CMR_MATRIX_INTERNAL_H

#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Abstract struct for row-wise representations of sparse matrices.
 */

typedef struct
{
  size_t numRows;       /**< \brief Number of rows. */
  size_t numColumns;    /**< \brief Number of columns. */
  size_t numNonzeros;   /**< \brief Number of and memory allocated for nonzeros. */
  size_t * rowSlice;    /**< \brief Array mapping each row to the index of its first entry. */
  size_t* entryColumns; /**< \brief Array mapping each entry to its column.*/
  void* entryValues;    /**< \brief Array mapping each entry to its value. */
} CMR_MATRIX;

/**
 * \brief Sorts the row and column indices of \p submatrix.
 */

CMR_ERROR CMRsortSubmatrix(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_SUBMAT* submatrix /**< The submatrix. */
);

/**
 * \brief Creates a \p numRows \f$ \times \f$ \p numColumns submatrix of the char \p matrix indexed by \p rows and
 * \p columns.
 */

CMR_ERROR CMRchrmatFilter(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Matrix. */
  size_t numRows,       /**< Number of rows of submatrix. */
  size_t* rows,         /**< Rows of \p matrix to be copied into the submatrix. */
  size_t numColumns,    /**< Number of columns of submatrix. */
  size_t* columns,      /**< Columns of \p matrix to be copied into the submatrix. */
  CMR_CHRMAT** presult  /**< Pointer for storing the created submatrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_MATRIX_INTERNAL_H */

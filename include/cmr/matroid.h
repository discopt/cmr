#ifndef CMR_MATROID_H
#define CMR_MATROID_H

/**
 * \file matroid.h
 *
 * \author Matthias Walter
 *
 * \brief Functionality for matroids.
 */

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup Matroid Matroid decomposition
 *
 * @{
 */

/**
 * \brief Applie a pivot to \p matrix and returns the resulting matrix in \p *presult.
 *
 * Calculations are done over the binary field.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatBinaryPivot(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Matrix to work with. */
  size_t pivotRow,      /**< Row of the pivot. */
  size_t pivotColumn,   /**< Column of the pivot. */
  CMR_CHRMAT** presult  /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief Applies a sequence of pivots to \p matrix and returns the resulting matrix in \p *presult.
 *
 * Calculations are done over the ternary field.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatTernaryPivot(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Matrix to work with. */
  size_t pivotRow,      /**< Row of the pivot. */
  size_t pivotColumn,   /**< Column of the pivot. */
  CMR_CHRMAT** presult  /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief Applies a sequence of pivots to \p matrix and returns the resulting matrix in \p *presult.
 *
 * Calculations are done over the binary field.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatBinaryPivots(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Matrix to work with. */
  size_t numPivots,       /**< Number of pivots to carry out. */
  size_t* pivotRows,      /**< Array with rows of the pivots. */
  size_t* pivotColumns,   /**< Array with columns of the pivots. */
  CMR_CHRMAT** presult    /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief Applies a sequence of pivots to \p matrix and returns the resulting matrix in \p *presult.
 *
 * Calculations are done over the ternary field.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatTernaryPivots(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Matrix to work with. */
  size_t numPivots,       /**< Number of pivots to carry out. */
  size_t* pivotRows,      /**< Array with rows of the pivots. */
  size_t* pivotColumns,   /**< Array with columns of the pivots. */
  CMR_CHRMAT** presult    /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief A minor of a matroid.
 *
 * Specified by a sequence of pivots and a submatrix.
 */

typedef struct
{
  size_t numPivots;               /**< Number of pivots to apply. */
  size_t* pivotRows;              /**< Array with pivot rows. */
  size_t* pivotColumns;           /**< Array with pivot columns. */
  CMR_SUBMAT* remainingSubmatrix; /**< Submatrix that one finally needs to look at. */
} CMR_MINOR;

/**
 * \brief Creates a minor, allocating space for \p numPivots pivots and a remaining \p submatrix.
 */

CMR_EXPORT
CMR_ERROR CMRminorCreate(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_MINOR** pminor,   /**< Pointer for storing the minor. */
  size_t numPivots,     /**< Number of pivots. */
  CMR_SUBMAT* submatrix /**< Submatrix (may be \c NULL; is not copied). */
);

/**
 * \brief Frees the minor \p *pminor (if \p pminor is not \c NULL).
 */

CMR_EXPORT
CMR_ERROR CMRminorFree(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_MINOR** pminor  /**< Pointer to the minor (may be \c NULL). */
);

/**
 * \brief Writes the minor \p minor to the file \fileName by means of lists of row and column indices as well as
 *        pivot entries.
 */

CMR_EXPORT
CMR_ERROR CMRminorPrint(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_MINOR* minor,   /**< Minor to write. */
  size_t numRows,     /**< Number of rows of original matrix. */
  size_t numColumns,  /**< Number of columns of original matrix. */
  FILE* stream        /**< File stream to save minor to.. */
);

/**
 * \brief Writes the minor \p minor to the file \fileName by means of lists of row and column indices as well as
 *        pivot entries.
 */

CMR_EXPORT
CMR_ERROR CMRminorWriteToFile(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_MINOR* minor,     /**< Minor to write. */
  size_t numRows,       /**< Number of rows of original matrix. */
  size_t numColumns,    /**< Number of columns of original matrix. */
  const char* fileName  /**< File name to save minor to; \c NULL indicates stdout. */
);


/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_MATROID_H */

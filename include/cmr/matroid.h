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
 * \brief Apply a pivot to \p matrix and returns the resulting matrix in \p *presult.
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

typedef enum
{
  CMR_MINOR_TYPE_DETERMINANT = -2,
    /**< A submatrix \f$ M \f$ with \f$ |\det(M)| = 2 \f$. */
  CMR_MINOR_TYPE_ENTRY = -1,
    /**< A submatrix \f$ M \f$ of size 1-by-1 exhibiting a bad entry. */
  CMR_MINOR_TYPE_CUSTOM = 0,
    /**< A custom minor. */
  CMR_MINOR_TYPE_U24 = 1,
    /**< A minor representing \f$ U^2_4 \f$. */
  CMR_MINOR_TYPE_FANO = 2,
    /**< A minor representing \f$ F_7 \f$. */
  CMR_MINOR_TYPE_FANO_DUAL = 3,
    /**< A minor representing \f$ F_7^\star \f$. */
  CMR_MINOR_TYPE_K5 = 4,
    /**< A minor representing \f$ M(K_5) \f$. */
  CMR_MINOR_TYPE_K5_DUAL = 5,
    /**< A minor representing \f$ M(K_5)^\star \f$. */
  CMR_MINOR_TYPE_K33 = 6,
    /**< A minor representing \f$ M(K_{3,3}) \f$. */
  CMR_MINOR_TYPE_K33_DUAL = 7,
    /**< A minor representing \f$ M(K_{3,3})^\star \f$. */
} CMR_MINOR_TYPE;

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
  CMR_MINOR_TYPE type;            /**< Type of minor. */
} CMR_MINOR;

/**
 * \brief Creates a minor, allocating space for \p numPivots pivots and a remaining \p submatrix.
 */

CMR_EXPORT
CMR_ERROR CMRminorCreate(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_MINOR** pminor,     /**< Pointer for storing the minor. */
  size_t numPivots,       /**< Number of pivots. */
  CMR_SUBMAT* submatrix,  /**< Submatrix (may be \c NULL; is not copied). */
  CMR_MINOR_TYPE type     /**< Type of minor. */
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
 * \brief Returns the type of \p minor.
 */

CMR_EXPORT
CMR_MINOR_TYPE CMRminorType(
  CMR_MINOR* minor  /**< Minor. */
);

/**
 * \brief Returns the number of pivots needed to make a \p minor visible.
 */

CMR_EXPORT
size_t CMRminorNumPivots(
  CMR_MINOR* minor  /**< Minor. */
);

/**
 * \brief Returns the array with pivot rows to make a \p minor visible.
 */

CMR_EXPORT
size_t* CMRminorPivotRows(
  CMR_MINOR* minor  /**< Minor. */
);

/**
 * \brief Returns the array with pivot columns to make a \p minor visible.
 */

CMR_EXPORT
size_t* CMRminorPivotColumns(
  CMR_MINOR* minor  /**< Minor. */
);

/**
 * \brief Returns the submatrix to take (after applying pivots) to make a \p minor visible.
 */

CMR_EXPORT
CMR_SUBMAT* CMRminorSubmatrix(
  CMR_MINOR* minor  /**< Minor. */
);

/**
 * \brief Writes the minor \p minor to \p stream by means of lists of row and column indices as well as pivot entries.
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
 * \brief Writes the minor \p minor to the file \p fileName by means of lists of row and column indices as well as
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

#ifndef CMR_SEPARATION_H
#define CMR_SEPARATION_H

#include <cmr/env.h>
#include <cmr/matrix.h>

/**
 * \file separation.h
 *
 * \author Matthias Walter
 *
 * \brief Data structures for k-separations.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  unsigned char* rowsToPart;        /**< \brief Indicates to which block each row belongs. Values above 1 are ignored. */
  unsigned char* columnsToPart;     /**< \brief Indicates to which block each column belongs. Values above 1 are ignored. */
  size_t numRows[2];                /**< \brief Indicates the number of rows of each part. */
  size_t numColumns[2];             /**< \brief Indicates the number of columns of each part. */
  size_t* rows[2];                  /**< \brief Array of sorted rows for each part. */
  size_t* columns[2];               /**< \brief Array of sorted columns for each part. */
  unsigned char rankRows0Columns1;  /**< \brief Rank of submatrix corresponding to rows in part 0 and columns in part 1. */
  unsigned char rankRows1Columns0;  /**< \brief Rank of submatrix corresponding to rows in part 1 and columns in part 0. */
  unsigned char* indicatorMemory;   /**< \brief Memory for \ref rowsToPart and \ref columnsToPart. */
  size_t* elementMemory;            /**< \brief Memory for \ref rows and \ref columns. */
} CMR_SEPA;

/**
 * \brief Creates a separation.
 *
 * Only the memory is allocated. The usualy way to initialize it is to fill the arrays \ref rowsToPart and
 * \ref columnsToPart and then call \ref CMRsepaInitialize.
 */

CMR_EXPORT
CMR_ERROR CMRsepaCreate(
  CMR* cmr,           /**< \ref CMR environment. */
  size_t numRows,     /**< Number of rows. */
  size_t numColumns,  /**< Number of columns. */
  CMR_SEPA** psepa    /**< Pointer for storing the created separation. */
);

/**
 * \brief Initializes a separation.
 *
 * Assumes that \p separation was created via \ref CMRsepaCreate and that all entries of \ref rowsToPart and
 * \ref columnsToPart are set to either 0 or 1.
 */

CMR_EXPORT
CMR_ERROR CMRsepaInitialize(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_SEPA* sepa,                   /**< Already created separation. */
  unsigned char rankRows0Columns1,  /**< Rank of submatrix corresponding to rows in part 0 and columns in part 1. */
  unsigned char rankRows1Columns0   /**< Rank of submatrix corresponding to rows in part 1 and columns in part 0. */
);

/**
 * \brief Frees a separation.
 */

CMR_EXPORT
CMR_ERROR CMRsepaFree(
  CMR* cmr,         /**< \ref CMR environment. */
  CMR_SEPA** psepa  /**< Pointer to separation. */
);

/**
 * \brief Checks for a given matrix whether the binary k-separation is also a ternary one.
 *
 * Checks, for a ternary input matrix \f$ M \f$ and a k-separation (\f$ k \in \{1,2,3\} \f$) of the (binary) support
 * matrix of \f$ M \f$, whether it is also a k-separation of \f$ M \f$ itself. The result is stored in \p *pisTernary.
 *
 * If the check fails, a certifying submatrix is returned.
 */

CMR_EXPORT
CMR_ERROR CMRsepaCheckTernary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEPA* sepa,         /**< Separation. */
  CMR_CHRMAT* matrix,     /**< Matrix. */
  bool* pisTernary,       /**< Pointer for storing whether the check passed. */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a violator submatrix (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_SEPARATION_H */

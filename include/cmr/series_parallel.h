#ifndef CMR_SERIES_PARALLEL_H
#define CMR_SERIES_PARALLEL_H

/**
 * \file series_parallel.h
 *
 * \author Matthias Walter
 *
 * \brief Recognition of [series-parallel matrices](\ref series-parallel).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/element.h>
#include <cmr/matrix.h>
#include <cmr/separation.h>

/**
 * \brief Statistics for series-parallel recognition algorithm.
 */

typedef struct
{
  uint32_t totalCount;      /**< Total number of invocations. */
  double totalTime;       /**< Total time of all invocations. */
  uint32_t reduceCount;     /**< Number of calls to reduction algorithm. */
  double reduceTime;      /**< Time of reduction algorithm calls. */
  uint32_t wheelCount;      /**< Number of wheel matrix searches. */
  double wheelTime;       /**< Time of wheel matrix searches. */
  uint32_t nonbinaryCount;  /**< Number of searches for \f$ M_2 \f$ matrix. */
  double nonbinaryTime;   /**< Time of searches for \f$ M_2 \f$ matrix. */
} CMR_SP_STATISTICS;

/**
 * \brief Initializes all statistics for series-parallel computations.
 */

CMR_EXPORT
CMR_ERROR CMRstatsSeriesParallelInit(
  CMR_SP_STATISTICS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for series-parallel computations.
 */

CMR_EXPORT
CMR_ERROR CMRstatsSeriesParallelPrint(
  FILE* stream,             /**< File stream to print to. */
  CMR_SP_STATISTICS* stats, /**< Pointer to statistics. */
  const char* prefix        /**< Prefix string to prepend to each printed line (may be \c NULL). */
);

/**
 * \brief Represents a series-parallel reduction
 */

typedef struct
{
  CMR_ELEMENT element; /**< Element (row/column) that is removed. */
  CMR_ELEMENT mate;    /**< Element is parallel to or in series with \ref element, or 0 for a zero row/column. */
} CMR_SP_REDUCTION;

/**
 * Prints the series-parallel \p reduction to \p buffer.
 */

CMR_EXPORT
char* CMRspReductionString(
  CMR_SP_REDUCTION reduction, /**< Series-parallel reduction. */
  char* buffer                /**< Buffer to write to.
                               **  If \c NULL, a static one is used which will be overwritten in the next call.
                               **  Otherwise, it must hold at least 51 bytes.
                               **/
);

/**
 * \brief Returns \c true if the series-parallel \p reduction removes a row, i.e., is series.
 */

static inline
bool CMRspIsRow(
  CMR_SP_REDUCTION reduction /**< Series-parallel reduction. */
)
{
  return reduction.element < 0;
}

/**
 * \brief Returns \c true if the series-parallel \p reduction removes a column, i.e., is parallel.
 */

static inline
bool CMRspIsColumn(
  CMR_SP_REDUCTION reduction /**< Series-parallel reduction. */
)
{
  return reduction.element > 0;
}

/**
 * \brief Returns \c true if the series-parallel \p reduction removes a zero vector.
 */

static inline
bool CMRspIsZero(
  CMR_SP_REDUCTION reduction /**< Series-parallel reduction. */
)
{
  return reduction.mate == 0;
}

/**
 * \brief Returns \c true if the series-parallel \p reduction removes a unit vector.
 */

static inline
bool CMRspIsUnit(
  CMR_SP_REDUCTION reduction /**< Series-parallel reduction. */
)
{
  return reduction.element * reduction.mate < 0;
}

/**
 * \brief Returns \c true if the series-parallel \p reduction removes a vector that is a copy of another vector.
 */

static inline
bool CMRspIsCopy(
  CMR_SP_REDUCTION reduction /**< Series-parallel reduction. */
)
{
  return reduction.element * reduction.mate > 0;
}

/**
 * \brief Returns \c true if the series-parallel \p reduction is valid.
 */

static inline
bool CMRspIsValid(
  CMR_SP_REDUCTION reduction /**< Series-parallel reduction. */
)
{
  return reduction.element != 0;
}

/**
 * \brief Finds all series-parallel reductions for the binary \p matrix \f$ A \f$.
 *
 * Let \f$ A \in \{-1,0,1\}^{m \times n} \f$ with \f$ k \f$ nonzeros.
 *
 * Denote by \f$ A' \f$ the (maximum) SP-reduced submatrix of \f$ A \f$.
 * If \p premainingSubmatrix is not \c NULL, then the rows and columns of \f$ A' \f$ are stored.
 *
 * If \p pviolatorSubmatrix is not \c NULL and \p matrix is not binary series-parallel, then a wheel-submatrix is
 * stored.
 *
 * The running time is \f$ \mathcal{O} (m + n + k) \f$ assuming no hashtable collisions.
 */

CMR_EXPORT
CMR_ERROR CMRtestBinarySeriesParallel(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Sparse char matrix. */
  bool* pisSeriesParallel,          /**< Pointer for storing the result. */
  CMR_SP_REDUCTION* reductions,     /**< Array for storing the SP-reductions. If not \c NULL, it must have
                                     **  capacity at least number of rows + number of columns. */
  size_t* pnumReductions,           /**< Pointer for storing the number of SP-reductions. */
  CMR_SUBMAT** preducedSubmatrix,   /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix,  /**< Pointer for storing a wheel-submatrix (may be \c NULL). */
  CMR_SP_STATISTICS* stats,         /**< Pointer to statistics (may be \c NULL). */
  double timeLimit                  /**< Time limit to impose. */
);

/**
 * \brief Finds all series-parallel reductions for the ternary \p matrix \f$ A \f$.
 *
 * Let \f$ A \in \{-1,0,+1\}^{m \times n} \f$ with \f$ k \f$ nonzeros.
 *
 * Denote by \f$ A' \f$ the (maximum) SP-reduced submatrix of \f$ A \f$.
 * If \p premainingSubmatrix is not \c NULL, then the rows and columns of \f$ A' \f$ are stored.
 *
 * If \p pviolatorSubmatrix is not \c NULL and \p matrix is not ternary series-parallel, then a signed wheel- or
 * \f$ M_2 \f$-submatrix is stored.
 *
 * The running time is \f$ \mathcal{O} (m + n + k) \f$ assuming no hashtable collisions.
 */

CMR_EXPORT
CMR_ERROR CMRtestTernarySeriesParallel(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Sparse char matrix. */
  bool* pisSeriesParallel,          /**< Pointer for storing the result. */
  CMR_SP_REDUCTION* reductions,     /**< Array for storing the SP-reductions. If not \c NULL, it must have
                                     **  capacity at least number of rows + number of columns. */
  size_t* pnumReductions,           /**< Pointer for storing the number of SP-reductions. */
  CMR_SUBMAT** preducedSubmatrix,   /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix,  /**< Pointer for storing a signed wheel- or \f$ M_2 \f$-submatrix (may be \c NULL). */
  CMR_SP_STATISTICS* stats,         /**< Pointer to statistics (may be \c NULL). */
  double timeLimit                  /**< Time limit to impose. */
);

/**
 * \brief Finds all series-parallel reductions for the binary \p matrix \f$ A \f$.
 *
 * Let \f$ A \in \{0,1\}^{m \times n} \f$ with \f$ k \f$ nonzeros.
 *
 * Denote by \f$ A' \f$ the (maximum) SP-reduced submatrix of \f$ A \f$.
 * If \p premainingSubmatrix is not \c NULL, then the rows and columns of \f$ A' \f$ are stored.
 *
 * If \p pviolatorSubmatrix is not \c NULL and \p matrix is not binary series-parallel, then a wheel-submatrix is
 * stored (unless a 2-separation is found; see below). Note that the row/column indices refer to \f$ A \f$.
 *
 * If \p pseparation is not \c NULL and during the search for a wheel-submatrix a 2-separation that does not correspond
 * to an SP reduction is found then such a 2-separation is returned and the algorithm terminates without returning a
 * wheel-submatrix. Note that \p *pseparation then contains row/column indices relative to \f$ A' \f$.
 *
 * The running time is \f$ \mathcal{O} (m + n + k) \f$ assuming no hashtable collisions.
 */

CMR_EXPORT
CMR_ERROR CMRdecomposeBinarySeriesParallel(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Sparse char matrix. */
  bool* pisSeriesParallel,          /**< Pointer for storing the result. */
  CMR_SP_REDUCTION* reductions,     /**< Array for storing the SP-reductions. If not \c NULL, it must have
                                     **  capacity at least number of rows + number of columns. */
  size_t maxNumReductions,          /**< Maximum number of SP-reductions. Stops when this would be exceeded. */
  size_t* pnumReductions,           /**< Pointer for storing the number of SP-reductions; stores \c SIZE_MAX if
                                     **< \p maxNumReductions was exceeded.  */
  CMR_SUBMAT** preducedSubmatrix,   /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix,  /**< Pointer for storing a wheel-submatrix (may be \c NULL). */
  CMR_SEPA** pseparation,           /**< Pointer for storing a 2-separation (may be \c NULL). */
  CMR_SP_STATISTICS* stats,         /**< Pointer to statistics (may be \c NULL). */
  double timeLimit                  /**< Time limit to impose. */
);


/**
 * \brief Finds all series-parallel reductions for the ternary \p matrix \f$ A \f$.
 *
 * Let \f$ A \in \{-1,0,+1\}^{m \times n} \f$ with \f$ k \f$ nonzeros.
 *
 * Denote by \f$ A' \f$ the (maximum) SP-reduced submatrix of \f$ A \f$.
 * If \p premainingSubmatrix is not \c NULL, then the rows and columns of \f$ A' \f$ are stored.
 *
 * If \p pviolatorSubmatrix is not \c NULL and \p matrix is not ternary series-parallel, then a signed wheel- or
 * \f$ M_2 \f$-submatrix is stored (unless a 2-separation is found; see below). Note that the row/column indices refer
 * to \f$ A \f$.
 *
 * If \p pseparation is not \c NULL and during the search for a signed wheel-submatrix a 2-separation that does not
 * correspond to an SP reduction is found then such a 2-separation is returned and the algorithm terminates without
 * returning a signed wheel- or \f$ M_2 \f$-submatrix. Note that \p *pseparation then contains row/column indices
 * relative to \f$ A' \f$.
 *
 * The running time is \f$ \mathcal{O} (m + n + k) \f$ assuming no hashtable collisions.
 *
 * \see CMRsubmatZoomSubmat() for turning \p *pviolatorSubmatrix into a submatrix of the SP-reduced submatrix.
 */

CMR_EXPORT
CMR_ERROR CMRdecomposeTernarySeriesParallel(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Sparse char matrix. */
  bool* pisSeriesParallel,          /**< Pointer for storing the result. */
  CMR_SP_REDUCTION* reductions,     /**< Array for storing the SP-reductions. If not \c NULL, it must have
                                     **  capacity at least number of rows + number of columns. */
  size_t maxNumReductions,          /**< Maximum number of SP-reductions. Stops when this would be exceeded. */
  size_t* pnumReductions,           /**< Pointer for storing the number of SP-reductions; stores \c SIZE_MAX if
                                     **< \p maxNumReductions was exceeded.  */
  CMR_SUBMAT** preducedSubmatrix,   /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix,  /**< Pointer for storing a signed wheel- or \f$ M_2 \f$-submatrix (may be
                                     **  \c NULL). */
  CMR_SEPA** pseparation,           /**< Pointer for storing a 2-separation (may be \c NULL). */
  CMR_SP_STATISTICS* stats,         /**< Pointer to statistics (may be \c NULL). */
  double timeLimit                  /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_SERIES_PARALLEL_H */

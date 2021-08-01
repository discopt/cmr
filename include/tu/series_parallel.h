#ifndef TU_SERIES_PARALLEL_H
#define TU_SERIES_PARALLEL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/element.h>
#include <tu/matrix.h>

/**
 * \brief Statistics for series-parallel recognition algorithm.
 */

typedef struct
{
  size_t totalCount;
  double totalTime;
  size_t reduceCount;
  double reduceTime;
  size_t wheelCount;
  double wheelTime;
} TU_SP_STATISTICS;

/**
 * \brief Represents a series-parallel operation
 */

typedef struct
{
  TU_ELEMENT element; /**< Element (row/column) that is removed. */
  TU_ELEMENT mate;    /**< Element is parallel to or in series with \ref element, or 0 for a zero row/column. */
} TU_SP_OPERATION;

/**
 * \brief Initializes all statistics for series-parallel computations.
 */

TU_EXPORT
TU_ERROR TUspInitStatistics(
  TU_SP_STATISTICS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for series-parallel computations.
 */

TU_EXPORT
TU_ERROR TUspPrintStatistics(
  FILE* stream,           /**< File stream to print to. */
  TU_SP_STATISTICS* stats /**< Pointer to statistics. */
);

/**
 * Prints the series-parallel \p operation to \p buffer.
 */

TU_EXPORT
char* TUspOperationString(
  TU_SP_OPERATION operation,  /**< Series-parallel operation. */
  char* buffer      /**< Buffer to write to.
                      * If \c NULL, a static one is used which will be overwritten in the next call.
                      * Otherwise, it must hold at least 51 bytes.
                      **/
);

/**
 * \brief Returns \c true if the series-parallel \p operation removes a row, i.e., is series.
 */

static inline
bool TUspIsRow(
  TU_SP_OPERATION operation /**< Series-parallel operation. */
)
{
  return operation.element < 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation removes a column, i.e., is parallel.
 */

static inline
bool TUspIsColumn(
  TU_SP_OPERATION operation /**< Series-parallel operation. */
)
{
  return operation.element > 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation removes a zero vector.
 */

static inline
bool TUspIsZero(
  TU_SP_OPERATION operation /**< Series-parallel operation. */
)
{
  return operation.mate == 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation removes a unit vector.
 */

static inline
bool TUspIsUnit(
  TU_SP_OPERATION operation /**< Series-parallel operation. */
)
{
  return operation.element * operation.mate < 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation removes a vector that is a copy of another vector.
 */

static inline
bool TUspIsCopy(
  TU_SP_OPERATION operation /**< Series-parallel operation. */
)
{
  return operation.element * operation.mate > 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation is valid.
 */

static inline
bool TUspIsValid(
  TU_SP_OPERATION operation /**< Series-parallel operation. */
)
{
  return operation.element != 0;
}

/**
 * \brief Finds all series-parallel reductions for the ternary \p matrix.
 *
 * If \p premainingSubmatrix is not \c NULL, then the SP-reduced submatrix is stored.
 *
 * If \p pwheelSubmatrix is not \c NULL and \p matrix is not series-parallel, then a wheel submatrix is stored.
 * This may cause overhead that is linear in the number of rows + number of columns + number
 *
 * If \p separationRank1Elements is not \c NULL, also \p pnumSeparationRank1Elements and \p pwheelSubmatrix must not be
 * \c NULL. If during the search for the wheel submatrix a 2-separation (not belonging to a series-parallel reduction)
 * is found, then this is stored and the algorithm terminates.
 *
 * If \p isSorted is \c true, then the running time is linear in the number of rows + number of columns + number of
 * nonzeros of \p matrix up to hashtable collisions.
 * Otherwise, extra overhead is caused by sorting all nonzeros.
 */

TU_EXPORT
TU_ERROR TUfindSeriesParallel(
  TU* tu,                               /**< \ref TU environment. */
  TU_CHRMAT* matrix,                    /**< Sparse char matrix. */
  TU_SP_OPERATION* operations,          /**< Array for storing the SP-reductions. Must be sufficiently large, e.g., number
                                         **< of rows + number of columns. */
  size_t* pnumOperations,               /**< Pointer for storing the number of SP-reductions. */
  TU_SUBMAT** preducedSubmatrix,        /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  TU_SUBMAT** pwheelSubmatrix,          /**< Pointer for storing a submatrix representing a wheel (may be \c NULL). */
  TU_ELEMENT* separationRank1Elements,  /**< Array for storing elements of one part of a 2-separation. If not \c NULL,
                                         **< it must have sufficient capacity. */
  size_t* pnumSeparationRank1Elements,  /**< Pointer for storing the number of elements stored in
                                         **< \p separationRank1Elements (may be \c NULL). */
  bool isSorted,                        /**< Whether the entries of \p matrix are sorted. */
  TU_SP_STATISTICS* stats               /**< Pointer to statistics (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SERIES_PARALLEL_H */

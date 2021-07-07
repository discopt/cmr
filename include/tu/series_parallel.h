#ifndef TU_SERIES_PARALLEL_H
#define TU_SERIES_PARALLEL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/element.h>
#include <tu/matrix.h>

/**
 * \brief Represents a series-parallel operation
 */  

typedef struct
{
  TU_ELEMENT element; /**< Element (row/column) that is removed. */
  TU_ELEMENT mate;    /**< Element is parallel to or in series with \ref element, or 0 for a zero row/column. */
} TU_SP;

/**
 * Prints the series-parallel \p operation to \p buffer.
 */

TU_EXPORT
char* TUspString(
  TU_SP operation,  /**< Series-parallel operation. */
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
  TU_SP operation /**< Series-parallel operation. */
)
{
  return operation.element < 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation removes a column, i.e., is parallel.
 */

static inline
bool TUspIsColumn(
    TU_SP operation /**< Series-parallel operation. */
)
{
  return operation.element > 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation removes a zero vector.
 */

static inline
bool TUspIsZero(
    TU_SP operation /**< Series-parallel operation. */
)
{
  return operation.mate == 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation removes a unit vector.
 */

static inline
bool TUspIsUnit(
    TU_SP operation /**< Series-parallel operation. */
)
{
  return operation.element * operation.mate < 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation removes a vector that is a copy of another vector.
 */

static inline
bool TUspIsCopy(
    TU_SP operation /**< Series-parallel operation. */
)
{
  return operation.element * operation.mate > 0;
}

/**
 * \brief Returns \c true if the series-parallel \p operation is valid.
 */

static inline
bool TUspIsValid(
    TU_SP operation /**< Series-parallel operation. */
)
{
  return operation.element != 0;
}

/**
 * \brief Finds all series or parallel elements of the ternary \p matrix.
 *
 * If \p isSorted is \c true, then the running time is linear in the number of rows + number of columns + number of
 * nonzeros of \p matrix.
 */

TU_EXPORT
TU_ERROR TUfindSeriesParallel(
  TU* tu,                           /**< \ref TU environment. */
  TU_CHRMAT* matrix,                /**< Sparse char matrix. */
  TU_SP* operations,                /**< Array for storing the operations. Must be sufficiently large. */
  size_t* pnumOperations,           /**< Pointer for storing the number of operations. */  
  TU_SUBMAT** premainingSubmatrix,  /**< Pointer for storing the submatrix that remains. */
//   TU_SUBMAT** pwheelSubmatrix,      /**< Pointer for storing a submatrix representing a wheel. */
  bool isSorted                     /**< Whether the entries of \p matrix are sorted. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SERIES_PARALLEL_H */

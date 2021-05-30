#ifndef TU_SERIES_PARALLEL_H
#define TU_SERIES_PARALLEL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/element.h>
#include <tu/matrix.h>

typedef struct
{
  TU_ELEMENT element;
  TU_ELEMENT mate;
} TU_SERIES_PARALLEL;

/**
 * \brief Splits off all series or parallel elements of \p matrix.
 *
 * The \p matrix is assumed to be ternary.
 */

TU_EXPORT
TU_ERROR TUsplitSeriesParallelChr(
  TU* tu,                         /**< \ref TU environment. */
  TU_CHRMAT* matrix,              /**< Sparse char matrix. */
  TU_SERIES_PARALLEL* operations, /**< Array for storing the operations. Must be sufficiently large. */
  size_t* pnumOperations          /**< Pointer for storing the number of operations. */  
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SERIES_PARALLEL_H */

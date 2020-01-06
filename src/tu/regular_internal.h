#ifndef TU_REGULAR_INTERNAL_H
#define TU_REGULAR_INTERNAL_H

#include <tu/regular.h>

/**
 * \brief Tests binary linear matroid for regularity. The given matrix must be sequentially
 *        connected.
 *
 * If \p decomposition is not \c NULL, \c *decomposition will be a decomposition tree.
 */

TU_EXPORT
bool TUtestBinaryRegularSequentiallyConnected(
  TU* tu,                     /**< TU environment */
  TU_MATRIX_CHAR* matrix,     /**< Double matrix */
  TU_MATRIX_CHAR* transpose,  /**< Double transpose */
  TU_DEC** decomposition,     /**< If not \c NULL, the decomposition tree is stored. */
  bool notGraphic,            /**< If \c true, the matroid is known to be not graphic. */
  bool notCographic           /**< If \c true, the matroid is known to be not cographic. */
);

#endif /* TU_REGULAR_INTERNAL_H */


#ifndef TU_REGULAR_INTERNAL_H
#define TU_REGULAR_INTERNAL_H

#include <tu/regular.h>

/**
 * \brief Creates a decomposition tree node.
 */
void TUcreateDec(
  TU* tu,       /**< TU environment. */
  TU_DEC** dec  /**< Pointer to decomposition tree. */
);

/**
 * \brief Tests binary linear matroid for regularity. The given decomposition is initialized with
 *        a sequentially connected matrix.
 */

bool TUregularSequentiallyConnected(
  TU* tu,                 /**< TU environment. */
  TU_DEC* decomposition,  /**< Partial decomposition. */
  bool certify,           /**< Whether to compute a proper decomposition tree. */
  bool notGraphic,        /**< If \c true, the matroid is known to be not graphic. */
  bool notCographic       /**< If \c true, the matroid is known to be not cographic. */
);

#endif /* TU_REGULAR_INTERNAL_H */


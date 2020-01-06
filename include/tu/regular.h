#ifndef TU_REGULAR_H
#define TU_REGULAR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/env.h>
#include <tu/matroid.h>
#include <tu/matrix.h>

/**
 * \brief Tests binary linear matroid for regularity.
 * 
 * If \p decomposition is not \c NULL, \c *decomposition will be a decomposition tree.
 */

TU_EXPORT
bool TUtestBinaryRegularLabeled(
  TU* tu,                 /**< TU environment */
  TU_CHAR_MATRIX* matrix, /**< Double matrix */
  int* rowLabels,         /**< Labels of matroid elements corresponding to rows. */
  int* columnLabels,      /**< Labels of matroid elements corresponding to columns. */
  TU_DEC** decomposition  /**< If not \c NULL, the decomposition tree is stored. */
);

/**
 * \brief Tests binary linear matroid for regularity.
 * 
 * If \p decomposition is not \c NULL, \c *decomposition will be a decomposition tree.
 */

TU_EXPORT
bool TUtestBinaryRegular(
  TU* tu,                 /**< TU environment */
  TU_CHAR_MATRIX* matrix, /**< Double matrix */
  TU_DEC** decomposition  /**< If not \c NULL, the decomposition tree is stored. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_REGULAR_H */


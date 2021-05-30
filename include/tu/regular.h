#ifndef TU_REGULAR_H
#define TU_REGULAR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/env.h>
#include <tu/matrix.h>
#include <tu/graph.h>
#include <tu/decomposition.h>

/**
 * \brief Tests binary linear matroid for regularity.
 *
 * If \p decomposition is not \c NULL, \c *decomposition will be a decomposition tree.
 *
 */

TU_EXPORT
TU_ERROR TUtestBinaryRegular(
  TU* tu,                   /**< \ref TU environment. */
  TU_CHRMAT* matrix,        /**< Char matrix. */
  TU_ELEMENT* rowElements,     /**< Row elements or \c NULL for canonical ones. */
  TU_ELEMENT* columnElements,  /**< Column elements or \c NULL for canonical ones. */
  bool checkPlanarity,      /**< Whether graphic minors should be checked for cographicness. */
  bool certify,             /**< Whether an \f$ F_7 \f$ or \f$ F_7^\star \F$ minor shall be identified. */
  bool *pisRegular,         /**< Whether \p matrix is regular. */
  TU_DEC** pdec             /**< Pointer for storing the decomposition tree (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_REGULAR_H */


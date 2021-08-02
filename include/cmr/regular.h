#ifndef CMR_REGULAR_H
#define CMR_REGULAR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/graph.h>
#include <cmr/decomposition.h>

/**
 * \brief Tests binary linear matroid for regularity.
 *
 * If \p decomposition is not \c NULL, \c *decomposition will be a decomposition tree.
 *
 */

CMR_EXPORT
CMR_ERROR CMRtestBinaryRegular(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,        /**< Char matrix. */
  CMR_ELEMENT* rowElements,     /**< Row elements or \c NULL for canonical ones. */
  CMR_ELEMENT* columnElements,  /**< Column elements or \c NULL for canonical ones. */
  bool checkPlanarity,      /**< Whether graphic minors should be checked for cographicness. */
  bool certify,             /**< Whether an \f$ F_7 \f$ or \f$ F_7^\star \F$ minor shall be identified. */
  bool *pisRegular,         /**< Whether \p matrix is regular. */
  CMR_TU_DEC** pdec             /**< Pointer for storing the decomposition tree (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_REGULAR_H */


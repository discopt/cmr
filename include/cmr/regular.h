#ifndef CMR_REGULAR_H
#define CMR_REGULAR_H

/**
 * \file regular.h
 *
 * \author Matthias Walter and Klaus Truemper
 *
 * \brief Recognition of [regular matrices](\ref regular).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/matroid.h>
#include <cmr/graph.h>
#include <cmr/dec.h>

/**
 * \brief Tests binary linear matroid for regularity.
 *
 * If \p pdec is not \c NULL, \c *pdec will be a (partial) decomposition tree.
 * If \p completeTree is \c true, then the decomposition tree is complete. Otherwise, it must only contain sufficient
 * information in order to determine regularity.
 *
 * If \p pminor is not \c NULL and \p matrix is not regular, then an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor is searched.
 * This causes additional computational effort!
 */

CMR_EXPORT
CMR_ERROR CMRtestBinaryRegular(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Input matrix. */
  bool *pisRegular,     /**< Pointer for storing whether \p matrix is regular. */
  CMR_DEC** pdec,       /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_MINOR** pminor,   /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  bool checkPlanarity,  /**< Whether graphic minors should be checked for cographicness. */
  bool completeTree     /**< Whether a complete decomposition tree shall be computed. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_REGULAR_H */


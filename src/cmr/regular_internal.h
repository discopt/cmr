#ifndef CMR_REGULAR_INTERNAL_H
#define CMR_REGULAR_INTERNAL_H

#include <cmr/regular.h>

/**
 * \brief Performs a 1-sum decomposition of \p matrix and stores it in \p dec.
 *
 * If \p matrix is 1-connected, then \p dec remains unchanged. In particular, the \c matrix and \c transpose members remain
 * \c NULL. Otherwise, \p dec will become a \ref CMR_TU_DEC_ONE_SUM node with children that are initialized to the
 * 1-connected components. In this case, the \c matrix and \c transpose members of the child nodes are set.
 */

CMR_ERROR CMRregularDecomposeOneSum(
  CMR* cmr,          /**< \ref CMR environment. */
  CMR_TU_DEC* dec,      /**< Initialized decomposition node corresponding to \p matrix. */
  CMR_CHRMAT* matrix /**< Matrix. */
);

/**
 * \brief Tests binary linear matroid for regularity. The given decomposition is initialized with a sequentially
 *        connected matrix.
 */

CMR_ERROR CMRtestBinaryRegularConnected(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_TU_DEC* dec,                 /**< Initialized decomposition node. */
  CMR_CHRMAT* matrix,           /**< Char matrix corresponding to \p dec. */
  CMR_CHRMAT* transpose,        /**< Transpose of \p matrix (may be \c NULL). */
  bool checkPlanarity,          /**< Whether graphic minors should be checked for cographicness. */
  bool certify,                 /**< Whether an \f$ F_7 \f$ or \f$ F_7^\star \F$ minor shall be identified. */
  bool constructDecomposition,  /**< Whether the construction of a decomposition is needed. */
  bool *pisRegular              /**< Whether \p matrix is regular. */
);

#endif /* CMR_REGULAR_INTERNAL_H */


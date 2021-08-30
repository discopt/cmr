#ifndef CMR_REGULAR_INTERNAL_H
#define CMR_REGULAR_INTERNAL_H

#include <cmr/regular.h>

/**
 * \brief Tests ternary or binary linear matroid for regularity.
 *
 * If \p pdec is not \c NULL, \c *pdec will be a (partial) decomposition tree.
 * If \p completeTree is \c true, then the decomposition tree is complete. Otherwise, it must only contain sufficient
 * information in order to determine regularity.
 *
 * If \p pminor is not \c NULL and \p matrix is not regular, then an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor or a
 * submatrix with non-ternary determinant is searched. This causes additional computational effort!
 */

CMR_ERROR CMRtestRegular(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Input matrix. */
  bool ternary,         /**< Whether signs of \p matrix play a role. */
  bool *pisRegular,     /**< Pointer for storing whether \p matrix is regular. */
  CMR_DEC** pdec,       /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_MINOR** pminor,   /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  bool checkPlanarity,  /**< Whether graphic minors should be checked for cographicness. */
  bool completeTree     /**< Whether a complete decomposition tree shall be computed. */
);

/**
 * \brief Performs a 1-sum decomposition of \p matrix and stores it in \p dec.
 *
 * If \p matrix is 1-connected, then \p dec remains unchanged. In particular, the \c matrix and \c transpose members remain
 * \c NULL. Otherwise, \p dec will become a \ref CMR_DEC_ONE_SUM node with children that are initialized to the
 * 1-connected components. In this case, the \c matrix and \c transpose members of the child nodes are set.
 */

CMR_ERROR CMRregularDecomposeOneSum(
  CMR* cmr,          /**< \ref CMR environment. */
  CMR_DEC* dec,      /**< Initialized decomposition node corresponding to \p matrix. */
  CMR_CHRMAT* matrix /**< Matrix. */
);

/**
 * \brief Splits off series-parallel elements from the matrix of a decomposition node.
 *
 * The decomposition node \p *pdec and matrix \p *pmatrix are replaced by the SP-reduced ones, i.e., by the
 * corresponding (grand-) children of the given decomposition and the corresponding matrix.
 */

CMR_ERROR CMRregularDecomposeSeriesParallel(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_DEC** pdec,           /**< Pointer to decomposition node. */
  CMR_CHRMAT** pmatrix,     /**< Pointer to matrix of \p *pdec. */
  bool ternary              /**< Whether to consider the signs of the matrix. */
);

#endif /* CMR_REGULAR_INTERNAL_H */

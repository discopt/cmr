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

typedef enum
{
  CMR_DEC_CONSTRUCT_NONE = 0,
  CMR_DEC_CONSTRUCT_LEAVES = 1,
  CMR_DEC_CONSTRUCT_ALL = 2,
} CMR_DEC_CONSTRUCT;

typedef struct
{
  bool planarityCheck;          /**< \brief Whether minors identified as graphic should still be checked for
                                 **         cographicness; default: \c false. */
  bool completeTree;            /**< \brief Whether to compute a complete decomposition tree (even if already
                                 **         non-regular; default: \c false. */
  CMR_DEC_CONSTRUCT matrices;   /**< \brief Which matrices of the decomposition to construct; default:
                                 **         \ref CMR_DEC_CONSTRUCT_NONE. */
  CMR_DEC_CONSTRUCT transposes; /**< \brief Which transposed matrices of the decomposition to construct; default:
                                 **         \ref CMR_DEC_CONSTRUCT_NONE. */
  CMR_DEC_CONSTRUCT graphs;     /**< \brief Which (co)graphs to construct; default: \ref CMR_DEC_CONSTRUCT_NONE. */
} CMR_REGULAR_PARAMETERS;

/**
 * \brief Initializes the default parameters for regularity testing.
 *
 * These are selected for minimum running time.
 */

CMR_EXPORT
CMR_ERROR CMRregularInitParameters(
  CMR_REGULAR_PARAMETERS* params  /**< Pointer to parameters. */
);
  
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
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Input matrix. */
  bool *pisRegular,               /**< Pointer for storing whether \p matrix is regular. */
  CMR_DEC** pdec,                 /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_MINOR** pminor,             /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  CMR_REGULAR_PARAMETERS* params  /**< Parameters for the computation (may be \c NULL for defaults). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_REGULAR_H */


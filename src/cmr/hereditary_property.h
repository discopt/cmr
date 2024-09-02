#ifndef CMR_HEREDITARY_PROPERTY_H
#define CMR_HEREDITARY_PROPERTY_H

#include <cmr/env.h>
#include "env_internal.h"
#include <cmr/matrix.h>
#include <cmr/element.h>

#ifdef __cplus
extern "C" {
#endif

typedef CMR_ERROR (*HereditaryPropertyTest)(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,       /**< Some matrix to be tested for the property. */
  void* data,               /**< Potential additional data for the test function. */
  bool* phasProperty,       /**< Pointer for storing whether \p matrix has the property. */
  CMR_SUBMAT** psubmatrix,  /**< Pointer for storing a proper submatrix of \p matrix without the property. */
  double timeLimit          /**< Time limit to impose. */
); /**< Function pointer for functions that test a hereditary matrix property. */

/**
 * \brief Tests a given \p matrix for the hereditary property defined by a given \p testFunction.
 *
 * The algorithm finds the submatrix by successively single zeroing out rows or columns.
 */

CMR_ERROR CMRtestHereditaryPropertyNaive(
  CMR* cmr,                             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,                   /**< Some matrix not having the hereditary property. */
  HereditaryPropertyTest testFunction,  /**< Test function. */
  void* testData,                       /**< Data to be forwarded to the test function. */
  CMR_SUBMAT** psubmatrix,              /**< Pointer for storing a minimal submatrix not having the property. */
  double timeLimit                      /**< Time limit to impose. */
);

/**
 * \brief Greedily tests a given \p matrix for the hereditary property defined by a given \p testFunction.
 *
 * The algorithm finds the submatrix by successively zeroing out sets of rows or columns.
 */

CMR_ERROR CMRtestHereditaryPropertyGreedy(
  CMR* cmr,                             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,                   /**< Some matrix not having the hereditary property. */
  HereditaryPropertyTest testFunction,  /**< Test function. */
  void* testData,                       /**< Data to be forwarded to the test function. */
  CMR_SUBMAT** psubmatrix,              /**< Pointer for storing a minimal submatrix not having the property. */
  double timeLimit                      /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_HEREDITARY_PROPERTY_H */

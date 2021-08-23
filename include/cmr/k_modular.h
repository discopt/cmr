#ifndef CMR_K_MODULAR_H
#define CMR_K_MODULAR_H

/**
 * \file k_modular.h
 *
 * \author Matthias Walter and Klaus Truemper
 *
 * \brief Recognition of [(strongly) k-modular matrices](\ref k-modular).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/matrix.h>

/**
 * \brief Tests a matrix \f$ M \f$ for being [unimodular](\ref k-modular).
 *
 * Tests if matrix \f$ M \f$ is unimodular and sets \p *pisUnimodular accordingly.
 */

CMR_EXPORT
CMR_ERROR CMRtestUnimodularity(
  CMR* cmr,           /**< \ref CMR environment */
  CMR_CHRMAT* matrix, /**< Matrix \f$ M \f$. */
  bool* pisUnimodular /**< Pointer for storing whether \f$ M \f$ is unimodular. */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [strongly unimodular](\ref k-modular).
 *
 * Tests if matrix \f$ M \f$ is strongly unimodular and sets \p *pisStronglyUnimodular accordingly.
 */

CMR_EXPORT
CMR_ERROR CMRtestStrongUnimodularity(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisStronglyUnimodular /**< Pointer for storing whether \f$ M \f$ is strongly unimodular. */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [k-modular](\ref k-modular).
 *
 * Tests if matrix \f$ M \f$ is k-modular and sets \p *pisKmodular accordingly.
 */

CMR_EXPORT
CMR_ERROR CMRtestKmodularity(
  CMR* cmr,           /**< \ref CMR environment */
  CMR_CHRMAT* matrix, /**< Matrix \f$ M \f$. */
  bool* pisKmodular,  /**< Pointer for storing whether \f$ M \f$ is k-modular. */
  size_t* pk          /**< Pointer for storing the value of \f$ k \f$. */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [strongly kmodular](\ref k-modular).
 *
 * Tests if matrix \f$ M \f$ is strongly k-modular and sets \p *pisStronglyKmodular accordingly.
 */

CMR_EXPORT
CMR_ERROR CMRtestStrongKmodularity(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisStronglyKmodular,  /**< Pointer for storing whether \f$ M \f$ is strongly k-modular. */
  size_t* pk                  /**< Pointer for storing the value of \f$ k \f$. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_UNIMODULAR_H */


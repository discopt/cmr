#ifndef CMR_CAMION_H
#define CMR_CAMION_H

/**
 * \file camion.h
 *
 * \author Matthias Walter
 *
 * \brief Testing whether a matrix is [Camion-signed](\ref camion).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/matrix.h>

/**
 * \brief Tests a matrix \f$ M \f$ for being a [Camion-signed](\ref camion).
 */

CMR_EXPORT
CMR_ERROR CMRtestCamionSigned(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Matrix \f$ M \f$. */
  bool* pisCamionSigned,  /**< Pointer for storing whether \f$ M \f$ is [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-Camion submatrix (may be \c NULL). */
);

/**
 * \brief Computes a [Camion-signed](\ref camion) version of a given ternary matrix \f$ M \f$.
 */

CMR_EXPORT
CMR_ERROR CMRcomputeCamionSigned(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Matrix \f$ M \f$ to be modified. */
  bool* pwasCamionSigned, /**< Pointer for storing whether \f$ M \f$ was already [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-Camion submatrix (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_CAMION_H */

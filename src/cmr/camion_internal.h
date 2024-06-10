#ifndef CMR_CAMION_INTERNAL_H
#define CMR_CAMION_INTERNAL_H

#include <cmr/env.h>
#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Ensures that sequentially connected matrix \f$ M \f$ is [Camion-signed](\ref camion).
 *
 * The matrix \f$ M \f$ is assumed to be ternary. If sign changes are necessary, only \p matrix is modified.
 * In particular, \p transpose remains unchanged.
 *
 * If \p submatrix is not \c NULL and sign changes are necessary, then a submatrix with determinant
 * -2 or +2 is stored in \p *psubmatrix and the caller must use \ref CMRsubmatFree() to free its
 * memory. It is set to \c NULL if no sign changes are needed.
 */

CMR_ERROR CMRcamionComputeSignSequentiallyConnected(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,       /**< Matrix \f$ M \f$. */
  CMR_CHRMAT* transpose,    /**< Transpose \f$ M^{\mathsf{T}} \f$. */
  bool change,              /**< Whether signs of \p matrix should be changed if necessary. */
  char* pmodification,      /**< Pointer for storing which matrix was modified. */
  CMR_SUBMAT** psubmatrix,  /**< Pointer for storing a submatrix with a bad determinant (may be \c NULL). */
  double timeLimit          /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_CAMION_INTERNAL_H */

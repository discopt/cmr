#ifndef CMR_INTERFACE_H
#define CMR_INTERFACE_H

#include <cmr/env.h>
#include <cmr/regular.h>
#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

CMR_ERROR CMRinterfaceTU(
  CMR* cmr,               /**< \ref CMR environment */
  CMR_CHRMAT* matrix,     /**< Matrix \f$ M \f$. */
  bool* pisTU,            /**< Pointer for storing whether \f$ M \f$ is TU. */
  CMR_TU_DEC** pdec,      /**< Pointer for storing the decomposition if \f$ M \f$ is TU (may be \c NULL). */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a bad submatrix with a bad determinant (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_INTERFACE_H */

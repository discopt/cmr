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
 * \brief Statistics for [Camion-signing](\ref camion) algorithm.
 */

typedef struct
{
  size_t totalCount;      /**< Total number of invocations. */
  double totalTime;       /**< Total time of all invocations. */
} CMR_CAMION_STATISTICS;

/**
 * \brief Initializes all statistics for [Camion-signing](\ref camion) algorithm.
 */

CMR_EXPORT
CMR_ERROR CMRstatsCamionInit(
  CMR_CAMION_STATISTICS* stats  /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for [Camion-signing](\ref camion) algorithm.
 */

CMR_EXPORT
CMR_ERROR CMRstatsCamionPrint(
  FILE* stream,                 /**< File stream to print to. */
  CMR_CAMION_STATISTICS* stats, /**< Pointer to statistics. */
  const char* prefix            /**< Prefix string to prepend to each printed line (may be \c NULL). */
);


/**
 * \brief Tests a matrix \f$ M \f$ for being a [Camion-signed](\ref camion).
 */

CMR_EXPORT
CMR_ERROR CMRtestCamionSigned(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Matrix \f$ M \f$. */
  bool* pisCamionSigned,        /**< Pointer for storing whether \f$ M \f$ is [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix,      /**< Pointer for storing a non-Camion submatrix (may be \c NULL). */
  CMR_CAMION_STATISTICS* stats  /**< Statistics for the computation (may be \c NULL). */
);

/**
 * \brief Computes a [Camion-signed](\ref camion) version of a given ternary matrix \f$ M \f$.
 */

CMR_EXPORT
CMR_ERROR CMRcomputeCamionSigned(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Matrix \f$ M \f$ to be modified. */
  bool* pwasCamionSigned,       /**< Pointer for storing whether \f$ M \f$ was already [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix,      /**< Pointer for storing a non-Camion submatrix (may be \c NULL). */
  CMR_CAMION_STATISTICS* stats  /**< Statistics for the computation (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_CAMION_H */

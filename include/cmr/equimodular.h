#ifndef CMR_EQUIMODULAR_H
#define CMR_EQUIMODULAR_H

/**
 * \file equimodular.h
 *
 * \author Matthias Walter and Klaus Truemper
 *
 * \brief Recognition of [(strongly) equimodular matrices](\ref equimodular).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/tu.h>

typedef struct
{
  CMR_TU_PARAMS tu; /**< \brief Parameters for TU test. */
} CMR_EQUIMODULAR_PARAMETERS;

/**
 * \brief Initializes the default parameters for recognition of [equimodular](\ref equimodular) matrices.
 *
 * These are selected for minimum running time.
 */

CMR_EXPORT
CMR_ERROR CMRparamsEquimodularityInit(
  CMR_EQUIMODULAR_PARAMETERS* params  /**< Pointer to parameters. */
);

/**
 * \brief Statistics for recognition algorithm for [equimodular](\ref equimodular) matrices.
 */

typedef struct
{
  uint32_t totalCount;  /**< Total number of invocations. */
  double totalTime;     /**< Total time of all invocations. */
  double linalgTime;    /**< Time spent on linear algebra. */
  CMR_TU_STATS tu;      /**< Total unimodularity test. */
} CMR_EQUIMODULAR_STATISTICS;

/**
 * \brief Initializes all statistics for recognition algorithm for [equimodular](\ref equimodular) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRstatsEquimodularityInit(
  CMR_EQUIMODULAR_STATISTICS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for recognition algorithm for [equimodular](\ref equimodular) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRstatsEquimodularityPrint(
  FILE* stream,                       /**< File stream to print to. */
  CMR_EQUIMODULAR_STATISTICS* stats,  /**< Pointer to statistics. */
  const char* prefix                  /**< Prefix string to prepend to each printed line (may be \c NULL). */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [equimodular](\ref equimodular) (for determinant gcd \f$ k \f$).
 *
 * Tests if matrix \f$ M \f$ is equimodular for determinant gcd \f$ k \f$ and sets \p *pisEquimodular accordingly.
 * If \p pgcdDet is not \c NULL, the behavior is as follows.
 * If \p *pgcdDet is positive, then it tests only for that particular value of \f$ k \f$.
 * Otherwise, \p *pgcdDet is set to \f$ k \f$ if \f$ M \f$ is equimodular for determinant gcd \f$ k \f$, and to \f$ 0 \f$ if \f$ M \f$ is not equimodular.
 */

CMR_EXPORT
CMR_ERROR CMRtestEquimodularity(
  CMR* cmr,                           /**< \ref CMR environment */
  CMR_INTMAT* matrix,                 /**< Matrix \f$ M \f$. */
  bool* pisEquimodular,               /**< Pointer for storing whether \f$ M \f$ is equimodular. */
  int64_t* pgcdDet,                   /**< Pointer for supplying/storing the determinant gcd (may be \c NULL). */
  CMR_EQUIMODULAR_PARAMETERS* params, /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_EQUIMODULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                    /**< Time limit to impose. */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [strongly equimodular](\ref equimodular).
 *
 * Tests if matrix \f$ M \f$ is strongly equimodular and sets \p *pisStronglyEquimodular accordingly.
 * If \p pgcdDet is not \c NULL, the behavior is as follows.
 * If \p *pgcdDet is positive, then it tests only for that particular value of \f$ k \f$.
 * Otherwise, \p *pgcdDet is set to \f$ k \f$ if \f$ M \f$ is strongly equimodular for determinant gcd \f$ k \f$,
 * and to \f$ 0 \f$ if \f$ M \f$ is not strongly equimodular.
 */

CMR_EXPORT
CMR_ERROR CMRtestStrongEquimodularity(
  CMR* cmr,                           /**< \ref CMR environment */
  CMR_INTMAT* matrix,                 /**< Matrix \f$ M \f$. */
  bool* pisStronglyEquimodular,       /**< Pointer for storing whether \f$ M \f$ is strongly equimodular. */
  int64_t* pgcdDet,                   /**< Pointer for supplying/storing the determinant gcd (may be \c NULL). */
  CMR_EQUIMODULAR_PARAMETERS* params, /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_EQUIMODULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                    /**< Time limit to impose. */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [unimodular](\ref equimodular).
 *
 * Tests if matrix \f$ M \f$ is unimodular and sets \p *pisUnimodular accordingly.
 */

CMR_EXPORT
CMR_ERROR CMRtestUnimodularity(
  CMR* cmr,                           /**< \ref CMR environment */
  CMR_INTMAT* matrix,                 /**< Matrix \f$ M \f$. */
  bool* pisUnimodular,                /**< Pointer for storing whether \f$ M \f$ is unimodular. */
  CMR_EQUIMODULAR_PARAMETERS* params, /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_EQUIMODULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                    /**< Time limit to impose. */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [strongly unimodular](\ref equimodular).
 *
 * Tests if matrix \f$ M \f$ is strongly unimodular and sets \p *pisStronglyUnimodular accordingly.
 */

CMR_EXPORT
CMR_ERROR CMRtestStrongUnimodularity(
  CMR* cmr,                           /**< \ref CMR environment */
  CMR_INTMAT* matrix,                 /**< Matrix \f$ M \f$. */
  bool* pisStronglyUnimodular,        /**< Pointer for storing whether \f$ M \f$ is strongly unimodular. */
  CMR_EQUIMODULAR_PARAMETERS* params, /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_EQUIMODULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                    /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_EQUIMODULAR_H */


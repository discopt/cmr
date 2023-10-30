#ifndef CMR_BALANCED_H
#define CMR_BALANCED_H

/**
 * \file balanced.h
 *
 * \author Henk Kraaij and Matthias Walter
 *
 * \brief Recognition of [balanced matrices](\ref balanced).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/camion.h>

/**
 * \brief Tests a matrix \f$ M \f$ for being [balanced](\ref balanced) via enumeration.
 *
 * Tests if matrix \f$ M \f$ is balanced via an enumeration algorithm and sets \p *pisBalanced accordingly.
 *
 * If \f$ M \f$ is not balanced and \p psubmatrix != \c NULL, then \p *psubmatrix will indicate a submatrix
 * of \f$ M \f$ with exactly two nonzeros in each row and in each column and with determinant \f$ -2 \f$ or \f$ 2 \f$.
 */

CMR_EXPORT
CMR_ERROR CMRtestBalancedEnumeration(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisBalanced,          /**< Pointer for storing whether \f$ M \f$ is balanced. */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a minimal nonbalanced submatrix (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [balanced](\ref balanced) via a graph-based algorithm.
 *
 * Tests if matrix \f$ M \f$ is balanced via a graph-based algorithm and sets \p *pisBalanced accordingly.
 *
 * If \f$ M \f$ is not balanced and \p psubmatrix != \c NULL, then \p *psubmatrix will indicate a submatrix
 * of \f$ M \f$ with exactly two nonzeros in each row and in each column and with determinant \f$ -2 \f$ or \f$ 2 \f$.
 */

CMR_EXPORT
CMR_ERROR CMRtestBalancedGraphBased(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisBalanced,          /**< Pointer for storing whether \f$ M \f$ is balanced. */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a minimal nonbalanced submatrix (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);
  
/**
 * \brief Tests a matrix \f$ M \f$ for being [balanced](\ref balanced).
 *
 * Tests if matrix \f$ M \f$ is balanced and sets \p *pisBalanced accordingly.
 * Automatically decides which algorithm to use.
 *
 * If \f$ M \f$ is not balanced and \p psubmatrix != \c NULL, then \p *psubmatrix will indicate a submatrix
 * of \f$ M \f$ with exactly two nonzeros in each row and in each column and with determinant \f$ -2 \f$ or \f$ 2 \f$.
 */

CMR_EXPORT
CMR_ERROR CMRtestBalanced(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisBalanced,          /**< Pointer for storing whether \f$ M \f$ is balanced. */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a minimal nonbalanced submatrix (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_BALANCED_H */

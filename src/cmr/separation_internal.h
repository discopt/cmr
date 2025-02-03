#ifndef CMR_SEPARATION_INTERNAL_H
#define CMR_SEPARATION_INTERNAL_H

#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Carries out the search for the connecting submatrix for a 3-sum.
 *
 * See \ref CMRthreesumDecomposeConnecting for the specification. This function carries out a breadth-first search in
 * the graphs of the top-left and/or bottom-right submatrices.
 *
 * At least one of \p topLeft and \p bottomRight must be \c true. If both are true then the implementation may detect a
 * submatrix with absolute determinant 2, which is also stored in \p *pviolator if \p pviolator is not \c NULL.
 * Even if this is the case, \p specialRows, \p specialColumns, \p *pgamma and \p *pbeta are set according to the
 * top-left matrix.
 *
 */

CMR_EXPORT
CMR_ERROR CMRthreesumDecomposeSearch(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Input matrix \f$ M \f$. */
  CMR_CHRMAT* transpose,  /**< Transpose matrix \f$ M^{\textsf{T}} \f$. */
  CMR_SEPA* sepa,         /**< 3-separation to decompose at. */
  bool topLeft,           /**< Whether to search in the top-left submatrix. */
  bool bottomRight,       /**< Whether to search in the bottom-right submatrix. */
  CMR_SUBMAT** pviolator, /**< Pointer for storing a submatrix with absolute determinant 2, if found. */
  size_t* specialRows,    /**< Array of length 2 for storing the rows \f$ i \f$ and \f$ j \f$ as rows of \f$ M \f$. */
  size_t* specialColumns, /**< Array of length 2 for storing the columns \f$ k \f$ and \f$ \ell \f$ as rows of
                           **  \f$ M \f$. */
  char* pgamma,           /**< Pointer for storing a correct value of \f$ \gamma \f$; may be \c NULL. */
  char* pbeta             /**< Pointer for storing a correct value of \f$ \beta \f$; may be \c NULL. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_SEPARATION_INTERNAL_H */

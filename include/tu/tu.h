#ifndef TU_TU_H
#define TU_TU_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/env.h>
#include <tu/matroid.h>
#include <tu/matrix.h>

/**
 * \brief Tests a double matrix for total unimodularity with absolute error
 * tolerance \p epsilon.
 * 
 * Returns \c true if and only if \p matrix is TU.
 *
 * If \p decomposition is not \c NULL and the algorithm has to test regularity of the support
 * matrix, then \c *decomposition will point to a decomposition tree for which the caller must use
 * \ref TUfreeDec to free memory. It is set to \c NULL otherwise.
 *
 * If \p submatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute
 * determinant larger than 1 will be searched, which may cause extra computational effort. In this
 * case, \c *submatrix will point to this submatrix for which the caller must use 
 * \ref TUfreeSubmatrix to free memory. It is set to \c NULL otherwise.
 */

TU_EXPORT
bool TUtestTotalUnimodularityDouble(
  TU* tu,                   /**< TU environment */
  TU_DOUBLE_MATRIX* matrix, /**< Double matrix */
  double epsilon,           /**< Absolute error tolerance */
  TU_DEC** decomposition,   /**< If not \c NULL, the decomposition tree is stored. */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

/**
 * \brief Tests an int matrix for total unimodularity.
 * 
 * Returns \c true if and only if \p matrix is TU.
 *
 * If \p decomposition is not \c NULL and the algorithm has to test regularity of the support
 * matrix, then \c *decomposition will point to a decomposition tree for which the caller must use
 * \ref TUfreeDec to free memory. It is set to \c NULL otherwise.
 *
 * If \p submatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute
 * determinant larger than 1 will be searched, which may cause extra computational effort. In this
 * case, \c *submatrix will point to this submatrix for which the caller must use 
 * \ref TUfreeSubmatrix to free memory. It is set to \c NULL otherwise.
 */

TU_EXPORT
bool TUtestTotalUnimodularityInt(
  TU* tu,                   /**< TU environment */
  TU_INT_MATRIX* matrix,    /**< Int matrix */
  TU_DEC** decomposition,   /**< If not \c NULL, the decomposition tree is stored. */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

/**
 * \brief Tests a char matrix for total unimodularity.
 * 
 * Returns \c true if and only if \p matrix is TU.
 *
 * If \p decomposition is not \c NULL and the algorithm has to test regularity of the support
 * matrix, then \c *decomposition will point to a decomposition tree for which the caller must use
 * \ref TUfreeDec to free memory. It is set to \c NULL otherwise.
 *
 * If \p submatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute
 * determinant larger than 1 will be searched, which may cause extra computational effort. In this
 * case, \c *submatrix will point to this submatrix for which the caller must use 
 * \ref TUfreeSubmatrix to free memory. It is set to \c NULL otherwise.
 */

TU_EXPORT
bool TUtestTotalUnimodularityChar(
  TU* tu,                   /**< TU environment */
  TU_CHAR_MATRIX* matrix,   /**< Char matrix */
  TU_DEC** decomposition,   /**< If not \c NULL, the decomposition tree is stored. */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_TU_H */

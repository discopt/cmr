#ifndef TU_TU_H
#define TU_TU_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/env.h>
#include <tu/matroid.h>

/**
 * \brief Tests if matrix is TU.
 * 
 * The dimensions of the matrix are \p numRows by \p numColumns. The arrays \p indexColumns and 
 * \p indexEntries of length \p numNonzeros contain column and value information, respectively, of
 * the nonzeros in row-wise order. The array \p rowStarts of length \p numRows contains the index of
 * the first entry of each row. The function returns \c true if and only if the matrix is TU.
 */

TU_EXPORT
bool TUtestTU(
  TU* tu, /**< TU environment */
  int numRows, /**< Number of matrix rows */
  int numColumns, /**< Number of matrix columns */
  int numNonzeros, /**< Number of matrix nonzeros */
  int* rowStarts, /**< Array with indices of first entries of each row. */
  int* indexColumns, /**< Array with columns for each index */
  char* indexEntries /**< Array with matrix entry for each index */
);

/**
 * \brief Tests if matrix is TU, returning a decomposition tree if applicable.
 * 
 * The dimensions of the matrix are \p numRows by \p numColumns. The arrays \p indexColumns and 
 * \p indexEntries of length \p numNonzeros contain column and value information, respectively, of
 * the nonzeros in row-wise order. The array \p rowStarts of length \p numRows contains the index of
 * the first entry of each row. The function returns \c true if and only if the matrix is TU.
 * 
 * If the algorithm has to test regularity of the support matrix, \p decomposition will point to a
 * decomposition tree, and is set to \c NULL otherwise.
 */

TU_EXPORT
bool TUtestTUwithDec(
  TU* tu, /**< TU environment */
  int numRows, /**< Number of matrix rows */
  int numColumns, /**< Number of matrix columns */
  int numNonzeros, /**< Number of matrix nonzeros */
  int* rowStarts, /**< Array with indices of first entries of each row. */
  int* indexColumns, /**< Array with columns for each index */
  char* indexEntries, /**< Array with matrix entry for each index */
  TU_DEC** decomposition /**< Decomposition tree */
);

/**
 * \brief Tests if matrix is TU, returning a violator if applicable.
 * 
 * The dimensions of the matrix are \p numRows by \p numColumns. The arrays \p indexColumns and 
 * \p indexEntries of length \p numNonzeros contain column and value information, respectively, of
 * the nonzeros in row-wise order. The array \p rowStarts of length \p numRows contains the index of
 * the first entry of each row. The function returns \c true if and only if the matrix is TU.
 * 
 * If the matrix is not TU, a violator, i.e., a submatrix with determinant -2 or 2 will be searched.
 * Note that this may cause extra computational effort. In this case, \p violator will point to this
 * submatrix, and is set to \c NULL otherwise.
 */

TU_EXPORT
bool TUtestTUwithSubmatrix(
  TU* tu, /**< TU environment */
  int numRows, /**< Number of matrix rows */
  int numColumns, /**< Number of matrix columns */
  int numNonzeros, /**< Number of matrix nonzeros */
  int* rowStarts, /**< Array with indices of first entries of each row. */
  int* indexColumns, /**< Array with columns for each index */
  char* indexEntries, /**< Array with matrix entry for each index */
  TU_SUBMATRIX** violator /**< Pointer to submatrix with determinant -2 or 2 */
);

/**
 * \brief Tests if matrix is TU, returning a decomposition tree or a violator if applicable.
 * 
 * The dimensions of the matrix are \p numRows by \p numColumns. The arrays \p indexColumns and 
 * \p indexEntries of length \p numNonzeros contain column and value information, respectively, of
 * the nonzeros in row-wise order. The array \p rowStarts of length \p numRows contains the index of
 * the first entry of each row. The function returns \c true if and only if the matrix is TU.
 * 
 * If \p decomposition is not \c NULL and if the algorithm has to test regularity of the support
 * matrix, then \p decomposition will point to a decomposition tree, and is set to \c NULL
 * otherwise.
 * 
 * If \p violator is not \c NULL and if the matrix is not TU, then a violator, i.e., a submatrix
 * with determinant -2 or 2 will be searched. Note that this may cause extra computational effort.
 * In this case, \p violator will point to this submatrix, and is set to \c NULL otherwise.
 */

TU_EXPORT
bool TUtestTUwithDecSubmatrix(
  TU* tu, /**< TU environment */
  int numRows, /**< Number of matrix rows */
  int numColumns, /**< Number of matrix columns */
  int numNonzeros, /**< Number of matrix nonzeros */
  int* rowStarts, /**< Array with indices of first entries of each row. */
  int* indexColumns, /**< Array with columns for each index */
  char* indexEntries, /**< Array with matrix entry for each index */
  TU_DEC** decomposition, /**< Decomposition tree */
  TU_SUBMATRIX** violator /**< Pointer to submatrix with determinant -2 or 2 */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_TU_H */

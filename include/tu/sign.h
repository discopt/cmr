#ifndef TU_SIGN_H
#define TU_SIGN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/env.h>

/**
 * \brief Tests if signs of matrix nonzeros qualify for being TU.
 *
 * The matrix is assumed to be ternary.
 *
 * The dimensions of the matrix are \p numRows by \p numColumns. The arrays \p indexColumns and 
 * \p indexEntries of length \p numNonzeros contain column and value information, respectively, of
 * the nonzeros in row-wise order. The array \p rowStarts of length \p numRows contains the index of
 * the first entry of each row.
 *
 * The function returns \c true if and only if the signs are correct.
 */

TU_EXPORT
bool TUtestSign(
  TU* tu, /**< TU environment */
  int numRows, /**< Number of matrix rows */
  int numColumns, /**< Number of matrix columns */
  int numNonzeros, /**< Number of matrix nonzeros */
  int* rowStarts, /**< Array with indices of first entries of each row. */
  int* indexColumns, /**< Array with columns for each index */
  char* indexEntries /**< Array with matrix entry for each index */
);


/**
 * \brief Modifies signs of matrix nonzeros to qualify for being TU.
 *
 * The matrix is assumed to be ternary.
 *
 * The dimensions of the matrix are \p numRows by \p numColumns. The arrays \p indexColumns and 
 * \p indexEntries of length \p numNonzeros contain column and value information, respectively, of
 * the nonzeros in row-wise order. The array \p rowStarts of length \p numRows contains the index of
 * the first entry of each row.
 *
 * The function returns \c true if and only if the signs were already correct.
 */

TU_EXPORT
bool TUcorrectSign(
  TU* tu, /**< TU environment */
  int numRows, /**< Number of matrix rows */
  int numColumns, /**< Number of matrix columns */
  int numNonzeros, /**< Number of matrix nonzeros */
  int* rowStarts, /**< Array with indices of first entries of each row. */
  int* indexColumns, /**< Array with columns for each index */
  char* indexEntries /**< Array with matrix entry for each index */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SIGN_H */

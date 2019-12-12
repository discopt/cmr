#include <tu/sign.h>

static bool sign(
  TU* tu, /**< TU environment */
  int numRows, /**< Number of matrix rows */
  int numColumns, /**< Number of matrix columns */
  int numNonzeros, /**< Number of matrix nonzeros */
  int* rowStarts, /**< Array with indices of first entries of each row. */
  int* indexColumns, /**< Array with columns for each index */
  char* indexEntries, /**< Array with matrix entry for each index */
  bool correct /**< Whether signs should be corrected */)
{
  
}

bool TUtestSign(TU* tu, int numRows, int numColumns, int numNonzeros, int* rowStarts,
  int* indexColumns, char* indexEntries)
{
  return sign(tu, numRows, numColumns, numNonzeros, rowStarts, indexColumns, indexEntries, false);
}

bool TUcorrectSign(TU* tu, int numRows, int numColumns, int numNonzeros, int* rowStarts,
  int* indexColumns, char* indexEntries)
{
  return sign(tu, numRows, numColumns, numNonzeros, rowStarts, indexColumns, indexEntries, true);
}

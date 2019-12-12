#include <tu/tu.h>

#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>

bool TUtestTU(
  TU* tu,
  int numRows,
  int numColumns,
  int numNonzeros,
  int* rowStarts,
  int* indexColumns,
  char* indexEntries
  )
{
  return TUtestTUwithDecSubmatrix(tu, numRows, numColumns, numNonzeros, rowStarts, indexColumns,
    indexEntries, NULL, NULL);
}

bool TUtestTUwithDec(
  TU* tu,
  int numRows,
  int numColumns,
  int numNonzeros,
  int* rowStarts,
  int* indexColumns,
  char* indexEntries,
  TU_DEC** decomposition
  )
{
  return TUtestTUwithDecSubmatrix(tu, numRows, numColumns, numNonzeros, rowStarts, indexColumns,
    indexEntries, decomposition, NULL);
}

bool TUtestTUwithSubmatrix(
  TU* tu,
  int numRows,
  int numColumns,
  int numNonzeros,
  int* rowStarts,
  int* indexColumns,
  char* indexEntries,
  TU_SUBMATRIX** violator
  )
{
  return TUtestTUwithDecSubmatrix(tu, numRows, numColumns, numNonzeros, rowStarts, indexColumns,
    indexEntries, NULL, violator);
}

bool TUtestTUwithDecSubmatrix(
  TU* tu,
  int numRows,
  int numColumns,
  int numNonzeros,
  int* rowStarts,
  int* indexColumns,
  char* indexEntries,
  TU_DEC** decomposition,
  TU_SUBMATRIX** violator
  )
{
  assert(false);
}

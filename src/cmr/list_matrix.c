#include "list_matrix.h"

#include <assert.h>
#include <stdint.h>

CMR_ERROR CMRlistmatrixCreate(CMR* cmr, size_t numRows, size_t numColumns, size_t memNonzeros, ListMatrix** presult)
{
  assert(cmr);
  assert(presult);

  CMR_CALL( CMRallocBlock(cmr, presult) );
  ListMatrix* result = *presult;

  result->numRows = numRows;
  result->rowElements = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &result->rowElements, numRows) );
  for (size_t row = 0; row < numRows; ++row)
  {
    result->rowElements[row].numNonzeros = 0;
    result->rowElements[row].head.row = row;
    result->rowElements[row].head.column = SIZE_MAX;
    result->rowElements[row].head.left = &result->rowElements[row].head;
    result->rowElements[row].head.right = &result->rowElements[row].head;
    result->rowElements[row].head.above = (row > 0) ? &result->rowElements[row - 1].head : &result->anchor;
    result->rowElements[row].head.below = (row + 1 < numRows) ? &result->rowElements[row + 1].head : &result->anchor;
    result->rowElements[row].head.value = 0;
    result->rowElements[row].head.special = 0;
  }

  result->numColumns = numColumns;
  result->columnElements = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &result->columnElements, numColumns) );
  for (size_t column = 0; column < numColumns; ++column)
  {
    result->columnElements[column].numNonzeros = 0;
    result->columnElements[column].head.column = column;
    result->columnElements[column].head.row = SIZE_MAX;
    result->columnElements[column].head.above = &result->columnElements[column].head;
    result->columnElements[column].head.below = &result->columnElements[column].head;
    result->columnElements[column].head.left =
      (column > 0) ? &result->columnElements[column - 1].head : &result->anchor;
    result->columnElements[column].head.right
      = (column + 1 < numColumns) ? &result->columnElements[column + 1].head : &result->anchor;
    result->columnElements[column].head.value = 0;
    result->columnElements[column].head.special = 0;
  }

  if (numRows > 0)
  {
    result->anchor.below = &result->rowElements[0].head;
    result->anchor.above = &result->rowElements[numRows - 1].head;
  }
  else
  {
    result->anchor.below = &result->anchor;
    result->anchor.above = &result->anchor;
  }

  if (numColumns > 0)
  {
    result->anchor.right = &result->columnElements[0].head;
    result->anchor.left = &result->columnElements[numColumns - 1].head;
  }
  else
  {
    result->anchor.right = &result->anchor;
    result->anchor.left = &result->anchor;
  }
  result->anchor.row = SIZE_MAX;
  result->anchor.column = SIZE_MAX;
  result->anchor.value = 0;
  result->anchor.special = 0;

  result->memNonzeros = memNonzeros;
  result->numNonzeros = 0;
  result->nonzeros = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &result->nonzeros, memNonzeros) );

  return CMR_OKAY;
}

CMR_ERROR CMRlistmatrixFree(CMR* cmr, ListMatrix ** plistmatrix)
{
  assert(cmr);
  assert(plistmatrix);

  ListMatrix* listmatrix = *plistmatrix;
  if (!listmatrix)
    return CMR_OKAY;

  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->nonzeros) );
  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->rowElements) );
  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->columnElements) );
  CMR_CALL( CMRfreeBlock(cmr, plistmatrix) );

  return CMR_OKAY;
}

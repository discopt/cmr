#include "listmatrix.h"

#include <assert.h>
#include <stdint.h>

CMR_ERROR CMRlistmatrixAlloc(CMR* cmr, size_t memRows, size_t memColumns, size_t memNonzeros, ListMatrix** presult)
{
  assert(cmr);
  assert(presult);

  CMR_CALL( CMRallocBlock(cmr, presult) );
  ListMatrix* result = *presult;

  result->numRows = 0;
  result->memRows = memRows;
  result->rowElements = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &result->rowElements, memRows) );

  result->numColumns = 0;
  result->memColumns = memColumns;
  result->columnElements = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &result->columnElements, memColumns) );
  
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

CMR_ERROR CMRlistmatrixInitializeFromMatrix(CMR* cmr, ListMatrix* listmatrix, CMR_CHRMAT* matrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);

  /* Reallocate if necessary. */
  if (listmatrix->memRows < matrix->numRows)
  {
    listmatrix->memRows = matrix->numRows;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->rowElements, matrix->numRows) );
  }
  if (listmatrix->memColumns < matrix->numColumns)
  {
    listmatrix->memColumns = matrix->numColumns;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->columnElements, matrix->numColumns) );
  }
  if (listmatrix->memNonzeros < matrix->numNonzeros)
  {
    listmatrix->memNonzeros = matrix->numNonzeros;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->nonzeros, matrix->numNonzeros) );
  }

  listmatrix->numRows = matrix->numRows;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    listmatrix->rowElements[row].numNonzeros = matrix->rowSlice[row + 1] - matrix->rowSlice[row];
    listmatrix->rowElements[row].head.row = row;
    listmatrix->rowElements[row].head.column = SIZE_MAX;
    listmatrix->rowElements[row].head.left = &listmatrix->rowElements[row].head;
    listmatrix->rowElements[row].head.right = &listmatrix->rowElements[row].head;
    listmatrix->rowElements[row].head.above = (row > 0) ? &listmatrix->rowElements[row - 1].head : &listmatrix->anchor;
    listmatrix->rowElements[row].head.below =
      (row + 1 < matrix->numRows) ? &listmatrix->rowElements[row + 1].head : &listmatrix->anchor;
    listmatrix->rowElements[row].head.value = 0;
    listmatrix->rowElements[row].head.special = 0;
  }

  listmatrix->numColumns = matrix->numColumns;
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    listmatrix->columnElements[column].numNonzeros = 0;
    listmatrix->columnElements[column].head.column = column;
    listmatrix->columnElements[column].head.row = SIZE_MAX;
    listmatrix->columnElements[column].head.above = &listmatrix->columnElements[column].head;
    listmatrix->columnElements[column].head.below = &listmatrix->columnElements[column].head;
    listmatrix->columnElements[column].head.left =
      (column > 0) ? &listmatrix->columnElements[column - 1].head : &listmatrix->anchor;
    listmatrix->columnElements[column].head.right
      = (column + 1 < matrix->numColumns) ? &listmatrix->columnElements[column + 1].head : &listmatrix->anchor;
    listmatrix->columnElements[column].head.value = 0;
    listmatrix->columnElements[column].head.special = 0;
  }

  if (matrix->numColumns > 0)
  {
    listmatrix->anchor.right = &listmatrix->columnElements[0].head;
    listmatrix->anchor.left = &listmatrix->columnElements[matrix->numColumns - 1].head;
  }
  else
  {
    listmatrix->anchor.right = &listmatrix->anchor;
    listmatrix->anchor.left = &listmatrix->anchor;
  }
  listmatrix->anchor.row = SIZE_MAX;
  listmatrix->anchor.column = SIZE_MAX;
  listmatrix->anchor.value = 0;
  listmatrix->anchor.special = 0;

  ListMatrixNonzero* nonzero = &listmatrix->nonzeros[0];
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      nonzero->row = row;
      nonzero->column = column;
      nonzero->value = matrix->entryValues[e];
      nonzero->special = 0;
      nonzero++;
      listmatrix->columnElements[column].numNonzeros++;
    }
  }

  /* Initialize linked list for rows. */
  if (matrix->numRows > 0)
  {
    listmatrix->anchor.below = &listmatrix->rowElements[0].head;
    listmatrix->rowElements[0].head.above = &listmatrix->anchor;
    listmatrix->anchor.above = &listmatrix->rowElements[matrix->numRows-1].head;
    listmatrix->rowElements[matrix->numRows-1].head.below = &listmatrix->anchor;
    for (size_t row = 1; row < matrix->numRows; ++row)
    {
      listmatrix->rowElements[row].head.above = &listmatrix->rowElements[row-1].head;
      listmatrix->rowElements[row-1].head.below = &listmatrix->rowElements[row].head;
    }
  }
  else
  {
    listmatrix->anchor.below = &listmatrix->anchor;
    listmatrix->anchor.above = &listmatrix->anchor;
  }

  /* Initialize linked list for columns. */
  if (matrix->numColumns > 0)
  {
    listmatrix->anchor.right = &listmatrix->columnElements[0].head;
    listmatrix->columnElements[0].head.left = &listmatrix->anchor;
    listmatrix->anchor.left = &listmatrix->columnElements[matrix->numColumns-1].head;
    listmatrix->columnElements[matrix->numColumns-1].head.right = &listmatrix->anchor;
    for (size_t column = 1; column < matrix->numColumns; ++column)
    {
      listmatrix->columnElements[column].head.left = &listmatrix->columnElements[column-1].head;
      listmatrix->columnElements[column-1].head.right = &listmatrix->columnElements[column].head;
    }
  }
  else
  {
    listmatrix->anchor.right = &listmatrix->anchor;
    listmatrix->anchor.left = &listmatrix->anchor;
  }

  /* Link the lists of nonzeros. */
  for (size_t i = 0; i < matrix->numNonzeros; ++i)
  {
    ListMatrixNonzero* nz = &listmatrix->nonzeros[i];

    nz->left = listmatrix->rowElements[listmatrix->nonzeros[i].row].head.left;
    listmatrix->rowElements[nz->row].head.left->right = nz;
    listmatrix->rowElements[nz->row].head.left = nz;

    nz->above = listmatrix->columnElements[listmatrix->nonzeros[i].column].head.above;
    listmatrix->columnElements[nz->column].head.above->below = nz;
    listmatrix->columnElements[nz->column].head.above = nz;
  }

  for (size_t row = 0; row < matrix->numRows; ++row)
    listmatrix->rowElements[row].head.left->right = &listmatrix->rowElements[row].head;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    listmatrix->columnElements[column].head.above->below = &listmatrix->columnElements[column].head;

  return CMR_OKAY;
}

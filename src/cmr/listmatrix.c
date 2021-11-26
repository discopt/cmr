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

CMR_ERROR CMRlistmatrixInitializeZero(CMR* cmr, ListMatrix* listmatrix, size_t numRows, size_t numColumns)
{
  assert(cmr);
  assert(listmatrix);

  /* Reallocate if necessary. */
  if (listmatrix->memRows < numRows)
  {
    listmatrix->memRows = numRows;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->rowElements, numRows) );
  }
  if (listmatrix->memColumns < numColumns)
  {
    listmatrix->memColumns = numColumns;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->columnElements, numColumns) );
  }

  /* Fill row data. */
  listmatrix->numRows = numRows;
  for (size_t row = 0; row < numRows; ++row)
  {
    listmatrix->rowElements[row].numNonzeros = 0;
    listmatrix->rowElements[row].head.row = row;
    listmatrix->rowElements[row].head.column = SIZE_MAX;
    listmatrix->rowElements[row].head.left = &listmatrix->rowElements[row].head;
    listmatrix->rowElements[row].head.right = &listmatrix->rowElements[row].head;
    listmatrix->rowElements[row].head.above = (row > 0) ? &listmatrix->rowElements[row - 1].head : &listmatrix->anchor;
    listmatrix->rowElements[row].head.below =
      (row + 1 < numRows) ? &listmatrix->rowElements[row + 1].head : &listmatrix->anchor;
    listmatrix->rowElements[row].head.value = 0;
    listmatrix->rowElements[row].head.special = 0;
  }

  /* Fill column data. */
  listmatrix->numColumns = numColumns;
  for (size_t column = 0; column < numColumns; ++column)
  {
    listmatrix->columnElements[column].numNonzeros = 0;
    listmatrix->columnElements[column].head.column = column;
    listmatrix->columnElements[column].head.row = SIZE_MAX;
    listmatrix->columnElements[column].head.above = &listmatrix->columnElements[column].head;
    listmatrix->columnElements[column].head.below = &listmatrix->columnElements[column].head;
    listmatrix->columnElements[column].head.left =
      (column > 0) ? &listmatrix->columnElements[column - 1].head : &listmatrix->anchor;
    listmatrix->columnElements[column].head.right
      = (column + 1 < numColumns) ? &listmatrix->columnElements[column + 1].head : &listmatrix->anchor;
    listmatrix->columnElements[column].head.value = 0;
    listmatrix->columnElements[column].head.special = 0;
  }

  /* Fill anchor data. */
  listmatrix->anchor.row = SIZE_MAX;
  listmatrix->anchor.column = SIZE_MAX;
  listmatrix->anchor.value = 0;
  listmatrix->anchor.special = 0;

  /* Initialize linked list for rows. */
  if (numRows > 0)
  {
    listmatrix->anchor.below = &listmatrix->rowElements[0].head;
    listmatrix->rowElements[0].head.above = &listmatrix->anchor;
    listmatrix->anchor.above = &listmatrix->rowElements[numRows-1].head;
    listmatrix->rowElements[numRows-1].head.below = &listmatrix->anchor;
    for (size_t row = 1; row < numRows; ++row)
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
  if (numColumns > 0)
  {
    listmatrix->anchor.right = &listmatrix->columnElements[0].head;
    listmatrix->columnElements[0].head.left = &listmatrix->anchor;
    listmatrix->anchor.left = &listmatrix->columnElements[numColumns-1].head;
    listmatrix->columnElements[numColumns-1].head.right = &listmatrix->anchor;
    for (size_t column = 1; column < numColumns; ++column)
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

  return CMR_OKAY;
}


CMR_ERROR CMRlistmatrixInitializeFromMatrix(CMR* cmr, ListMatrix* listmatrix, CMR_CHRMAT* matrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);

  /* Reallocate if necessary. */
  if (listmatrix->memNonzeros < matrix->numNonzeros)
  {
    listmatrix->memNonzeros = matrix->numNonzeros;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->nonzeros, matrix->numNonzeros) );
  }
  listmatrix->numNonzeros = matrix->numNonzeros;

  /* Initialze the zero matrix. */
  CMR_CALL( CMRlistmatrixInitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  /* Fill nonzero data. */
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
      listmatrix->rowElements[row].numNonzeros++;
      listmatrix->columnElements[column].numNonzeros++;
    }
  }

  /* Link the lists of nonzeros (left and above pointers). */
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

  /* Set the right and below pointers for nonzeros. */
  for (size_t row = 0; row < matrix->numRows; ++row)
    listmatrix->rowElements[row].head.left->right = &listmatrix->rowElements[row].head;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    listmatrix->columnElements[column].head.above->below = &listmatrix->columnElements[column].head;

  return CMR_OKAY;
}

CMR_ERROR CMRlistmatrixInitializeFromSubmatrix(CMR* cmr, ListMatrix* listmatrix, CMR_CHRMAT* matrix,
  CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(submatrix);

  CMR_CALL( CMRlistmatrixInitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  bool* rowUsed = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowUsed, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    rowUsed[row] = false;
  for (size_t r = 0; r < submatrix->numRows; ++r)
    rowUsed[submatrix->rows[r]] = true;

  bool* columnUsed = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnUsed, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnUsed[column] = false;
  for (size_t c = 0; c < submatrix->numColumns; ++c)
    columnUsed[submatrix->columns[c]] = true;

  /* Fill nonzero data. */
  ListMatrixNonzero* nonzero = &listmatrix->nonzeros[0];
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (!rowUsed[row])
      continue;

    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      if (!columnUsed[column])
        continue;

      nonzero->row = row;
      nonzero->column = column;
      nonzero->value = matrix->entryValues[e];
      nonzero->special = 0;
      listmatrix->numNonzeros++;

      nonzero->left = listmatrix->rowElements[row].head.left;
      listmatrix->rowElements[row].head.left->right = nonzero;
      nonzero->right = &listmatrix->rowElements[row].head;
      listmatrix->rowElements[row].head.left = nonzero;
      listmatrix->rowElements[row].numNonzeros++;

      nonzero->above = listmatrix->columnElements[column].head.above;
      listmatrix->columnElements[column].head.above->below = nonzero;
      nonzero->below = &listmatrix->columnElements[column].head;
      listmatrix->columnElements[column].head.above = nonzero;
      listmatrix->columnElements[column].numNonzeros++;

      nonzero++;
    }
  }

  CMR_CALL( CMRfreeStackArray(cmr, &columnUsed) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowUsed) );


  return CMR_OKAY;
}

CMR_ERROR CMRlistmatrixInitializeFromSubmatrixComplement(CMR* cmr, ListMatrix* listmatrix, CMR_CHRMAT* matrix,
  CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(submatrix);

  CMR_CALL( CMRlistmatrixInitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  bool* rowUsed = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowUsed, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    rowUsed[row] = false;
  for (size_t r = 0; r < submatrix->numRows; ++r)
    rowUsed[submatrix->rows[r]] = true;

  bool* columnUsed = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnUsed, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnUsed[column] = false;
  for (size_t c = 0; c < submatrix->numColumns; ++c)
    columnUsed[submatrix->columns[c]] = true;

  /* Fill nonzero data. */
  ListMatrixNonzero* nonzero = &listmatrix->nonzeros[0];
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      if (rowUsed[row] && columnUsed[column])
        continue;

      nonzero->row = row;
      nonzero->column = column;
      nonzero->value = matrix->entryValues[e];
      nonzero->special = 0;
      listmatrix->numNonzeros++;

      nonzero->left = listmatrix->rowElements[row].head.left;
      listmatrix->rowElements[row].head.left->right = nonzero;
      nonzero->right = &listmatrix->rowElements[row].head;
      listmatrix->rowElements[row].head.left = nonzero;
      listmatrix->rowElements[row].numNonzeros++;

      nonzero->above = listmatrix->columnElements[column].head.above;
      listmatrix->columnElements[column].head.above->below = nonzero;
      nonzero->below = &listmatrix->columnElements[column].head;
      listmatrix->columnElements[column].head.above = nonzero;
      listmatrix->columnElements[column].numNonzeros++;

      nonzero++;
    }
  }

  CMR_CALL( CMRfreeStackArray(cmr, &columnUsed) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowUsed) );

  return CMR_OKAY;
}

CMR_ERROR CMRlistmatrixPrintDense(CMR* cmr, ListMatrix* listmatrix, FILE* stream)
{
  assert(cmr);
  assert(listmatrix);
  assert(stream);

  for (size_t row = 0; row < listmatrix->numRows; ++row)
  {
    ListMatrixNonzero* nz = listmatrix->rowElements[row].head.right;
    for (size_t column = 0; column < listmatrix->numColumns; ++column)
    {
      if (column < nz->column)
        fprintf(stream, " 0");
      else
      {
        assert(column == nz->column);
        fprintf(stream, " %d", nz->value);
        nz = nz->right;
      }
      fflush(stream);
    }
    fprintf(stream, "\n");
  }
  fflush(stream);

  return CMR_OKAY;
}

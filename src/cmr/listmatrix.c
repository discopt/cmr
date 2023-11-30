// #define CMR_DEBUG

#include "listmatrix.h"

#include <assert.h>
#include <stdint.h>
#include <math.h>

CMR_ERROR CMRlistmat8Alloc(CMR* cmr, size_t memRows, size_t memColumns, size_t memNonzeros, ListMat8** presult)
{
  assert(cmr);
  assert(presult);

  CMR_CALL( CMRallocBlock(cmr, presult) );
  ListMat8* result = *presult;

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

CMR_ERROR CMRlistmat64Alloc(CMR* cmr, size_t memRows, size_t memColumns, size_t memNonzeros, ListMat64** presult)
{
  assert(cmr);
  assert(presult);

  CMR_CALL( CMRallocBlock(cmr, presult) );
  ListMat64* result = *presult;

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

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPAlloc(CMR* cmr, size_t memRows, size_t memColumns, size_t memNonzeros, ListMatGMP** presult)
{
  assert(cmr);
  assert(presult);

  CMR_CALL( CMRallocBlock(cmr, presult) );
  ListMatGMP* result = *presult;

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
  for (size_t i = 0; i < memNonzeros; ++i)
    mpz_init(result->nonzeros[i].value);

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8Free(CMR* cmr, ListMat8** plistmatrix)
{
  assert(cmr);
  assert(plistmatrix);

  ListMat8* listmatrix = *plistmatrix;
  if (!listmatrix)
    return CMR_OKAY;

  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->nonzeros) );
  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->rowElements) );
  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->columnElements) );
  CMR_CALL( CMRfreeBlock(cmr, plistmatrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64Free(CMR* cmr, ListMat64** plistmatrix)
{
  assert(cmr);
  assert(plistmatrix);

  ListMat64* listmatrix = *plistmatrix;
  if (!listmatrix)
    return CMR_OKAY;

  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->nonzeros) );
  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->rowElements) );
  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->columnElements) );
  CMR_CALL( CMRfreeBlock(cmr, plistmatrix) );

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPFree(CMR* cmr, ListMatGMP** plistmatrix)
{
  assert(cmr);
  assert(plistmatrix);

  ListMatGMP* listmatrix = *plistmatrix;
  if (!listmatrix)
    return CMR_OKAY;

  for (size_t i = 0; i < listmatrix->memNonzeros; ++i)
    mpz_clear(listmatrix->nonzeros[i].value);

  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->nonzeros) );
  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->rowElements) );
  CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->columnElements) );
  CMR_CALL( CMRfreeBlock(cmr, plistmatrix) );

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8InitializeZero(CMR* cmr, ListMat8* listmatrix, size_t numRows, size_t numColumns)
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

  /* Initialize the free list. */
  if (listmatrix->memNonzeros > 0)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[0];
    for (size_t i = 0; i < listmatrix->memNonzeros - 1; ++i)
    {
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
      listmatrix->nonzeros[i].value = 0;
    }
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
    listmatrix->nonzeros[listmatrix->memNonzeros-1].value = 0;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64InitializeZero(CMR* cmr, ListMat64* listmatrix, size_t numRows, size_t numColumns)
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

  /* Initialize the free list. */
  if (listmatrix->memNonzeros > 0)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[0];
    for (size_t i = 0; i < listmatrix->memNonzeros - 1; ++i)
    {
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
      listmatrix->nonzeros[i].value = 0;
    }
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
    listmatrix->nonzeros[listmatrix->memNonzeros-1].value = 0;
  }

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPInitializeZero(CMR* cmr, ListMatGMP* listmatrix, size_t numRows, size_t numColumns)
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
    mpz_init(listmatrix->rowElements[row].head.value);
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
    mpz_init(listmatrix->columnElements[column].head.value);
    listmatrix->columnElements[column].head.special = 0;
  }

  /* Fill anchor data. */
  listmatrix->anchor.row = SIZE_MAX;
  listmatrix->anchor.column = SIZE_MAX;
  mpz_init(listmatrix->anchor.value);
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

  /* Initialize the free list. */
  if (listmatrix->memNonzeros > 0)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[0];
    for (size_t i = 0; i < listmatrix->memNonzeros - 1; ++i)
    {
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
      mpz_init(listmatrix->nonzeros[i].value);
    }
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
    mpz_init(listmatrix->nonzeros[listmatrix->memNonzeros-1].value);
  }

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8InitializeFromChrMatrix(CMR* cmr, ListMat8* listmatrix, CMR_CHRMAT* matrix)
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
  CMR_CALL( CMRlistmat8InitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  /* Fill nonzero data. */
    ListMat8Nonzero* nonzero = &listmatrix->nonzeros[0];
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
    ListMat8Nonzero* nz = &listmatrix->nonzeros[i];

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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64InitializeFromIntMatrix(CMR* cmr, ListMat64* listmatrix, CMR_INTMAT* matrix)
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
  CMR_CALL( CMRlistmat64InitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  /* Fill nonzero data. */
  ListMat64Nonzero* nonzero = &listmatrix->nonzeros[0];
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
    ListMat64Nonzero* nz = &listmatrix->nonzeros[i];

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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPInitializeFromIntMatrix(CMR* cmr, ListMatGMP* listmatrix, CMR_INTMAT* matrix)
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
  CMR_CALL( CMRlistmatGMPInitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  /* Fill nonzero data. */
  ListMatGMPNonzero* nonzero = &listmatrix->nonzeros[0];
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      nonzero->row = row;
      nonzero->column = column;
      mpz_set_si(nonzero->value, matrix->entryValues[e]);
      nonzero->special = 0;
      nonzero++;
      listmatrix->rowElements[row].numNonzeros++;
      listmatrix->columnElements[column].numNonzeros++;
    }
  }

  /* Link the lists of nonzeros (left and above pointers). */
  for (size_t i = 0; i < matrix->numNonzeros; ++i)
  {
    ListMatGMPNonzero* nz = &listmatrix->nonzeros[i];

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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8InitializeFromDoubleMatrix(CMR* cmr, ListMat8* listmatrix, CMR_DBLMAT* matrix, double epsilon)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(epsilon >= 0);

  /* Reallocate if necessary. */
  if (listmatrix->memNonzeros < matrix->numNonzeros)
  {
    listmatrix->memNonzeros = matrix->numNonzeros;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->nonzeros, matrix->numNonzeros) );
  }
  listmatrix->numNonzeros = matrix->numNonzeros;

  /* Initialze the zero matrix. */
  CMR_CALL( CMRlistmat8InitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  /* Fill nonzero data. */
    ListMat8Nonzero* nonzero = &listmatrix->nonzeros[0];
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      nonzero->row = row;
      nonzero->column = column;
      double rounded = round(matrix->entryValues[e]);
      if (rounded > 127 || rounded < -127 || fabs(rounded - matrix->entryValues[e]) > epsilon)
        nonzero->value = -128;
      else
      {
        nonzero->value = (int8_t) rounded;
        assert(nonzero->value != -128);
      }
      nonzero->special = 0;
      nonzero++;
      listmatrix->rowElements[row].numNonzeros++;
      listmatrix->columnElements[column].numNonzeros++;
    }
  }

  /* Link the lists of nonzeros (left and above pointers). */
  for (size_t i = 0; i < matrix->numNonzeros; ++i)
  {
        ListMat8Nonzero* nz = &listmatrix->nonzeros[i];

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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64InitializeFromDoubleMatrix(CMR* cmr, ListMat64* listmatrix, CMR_DBLMAT* matrix, double epsilon)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(epsilon >= 0);

  /* Reallocate if necessary. */
  if (listmatrix->memNonzeros < matrix->numNonzeros)
  {
    listmatrix->memNonzeros = matrix->numNonzeros;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->nonzeros, matrix->numNonzeros) );
  }
  listmatrix->numNonzeros = matrix->numNonzeros;

  /* Initialze the zero matrix. */
  CMR_CALL( CMRlistmat64InitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  /* Fill nonzero data. */
  ListMat64Nonzero* nonzero = &listmatrix->nonzeros[0];
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      nonzero->row = row;
      nonzero->column = column;
      double rounded = round(matrix->entryValues[e]);
      if (rounded > INT32_MAX || rounded < -INT32_MAX || fabs(rounded - matrix->entryValues[e]) > epsilon)
        nonzero->value = INT32_MIN;
      else
      {
        nonzero->value = (int8_t) rounded;
        assert(nonzero->value != -128);
      }
      nonzero->special = 0;
      nonzero++;
      listmatrix->rowElements[row].numNonzeros++;
      listmatrix->columnElements[column].numNonzeros++;
    }
  }

  /* Link the lists of nonzeros (left and above pointers). */
  for (size_t i = 0; i < matrix->numNonzeros; ++i)
  {
    ListMat64Nonzero* nz = &listmatrix->nonzeros[i];

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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPInitializeFromDoubleMatrix(CMR* cmr, ListMatGMP* listmatrix, CMR_DBLMAT* matrix, double epsilon)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(epsilon >= 0);

  /* Reallocate if necessary. */
  if (listmatrix->memNonzeros < matrix->numNonzeros)
  {
    listmatrix->memNonzeros = matrix->numNonzeros;
    CMR_CALL( CMRreallocBlockArray(cmr, &listmatrix->nonzeros, matrix->numNonzeros) );
  }
  listmatrix->numNonzeros = matrix->numNonzeros;

  /* Initialze the zero matrix. */
  CMR_CALL( CMRlistmatGMPInitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

  /* Fill nonzero data. */
  ListMatGMPNonzero* nonzero = &listmatrix->nonzeros[0];
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      nonzero->row = row;
      nonzero->column = column;
      double rounded = round(matrix->entryValues[e]);
      if (rounded > 127 || rounded < -127 || fabs(rounded - matrix->entryValues[e]) > epsilon)
        mpz_init_set_si(nonzero->value, -128);
      else
      {
        mpz_init_set_si(nonzero->value, (int8_t) rounded);
        assert(mpz_cmp_si(nonzero->value, -128) != 0);
      }
      nonzero->special = 0;
      nonzero++;
      listmatrix->rowElements[row].numNonzeros++;
      listmatrix->columnElements[column].numNonzeros++;
    }
  }

  /* Link the lists of nonzeros (left and above pointers). */
  for (size_t i = 0; i < matrix->numNonzeros; ++i)
  {
    ListMatGMPNonzero* nz = &listmatrix->nonzeros[i];

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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8InitializeFromChrSubmatrix(CMR* cmr, ListMat8* listmatrix, CMR_CHRMAT* matrix,
  CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(submatrix);

  CMR_CALL( CMRlistmat8InitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

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
  ListMat8Nonzero* nonzero = &listmatrix->nonzeros[0];
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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros; ++i)
    {
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    }
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64InitializeFromIntSubmatrix(CMR* cmr, ListMat64* listmatrix, CMR_INTMAT* matrix,
  CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(submatrix);

  CMR_CALL( CMRlistmat64InitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

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
  ListMat64Nonzero* nonzero = &listmatrix->nonzeros[0];
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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPInitializeFromIntSubmatrix(CMR* cmr, ListMatGMP* listmatrix, CMR_INTMAT* matrix,
  CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(submatrix);

  CMR_CALL( CMRlistmatGMPInitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

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
  ListMatGMPNonzero* nonzero = &listmatrix->nonzeros[0];
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
      mpz_init_set_si(nonzero->value, matrix->entryValues[e]);
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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8InitializeFromChrSubmatrixComplement(CMR* cmr, ListMat8* listmatrix, CMR_CHRMAT* matrix,
  CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(submatrix);

  CMR_CALL( CMRlistmat8InitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

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
  ListMat8Nonzero* nonzero = &listmatrix->nonzeros[0];
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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64InitializeFromIntSubmatrixComplement(CMR* cmr, ListMat64* listmatrix, CMR_INTMAT* matrix,
  CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(submatrix);

  CMR_CALL( CMRlistmat64InitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

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
  ListMat64Nonzero* nonzero = &listmatrix->nonzeros[0];
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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPInitializeFromIntSubmatrixComplement(CMR* cmr, ListMatGMP* listmatrix, CMR_INTMAT* matrix,
  CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(listmatrix);
  assert(matrix);
  assert(submatrix);

  CMR_CALL( CMRlistmatGMPInitializeZero(cmr, listmatrix, matrix->numRows, matrix->numColumns) );

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
  ListMatGMPNonzero* nonzero = &listmatrix->nonzeros[0];
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
      mpz_init_set_si(nonzero->value, matrix->entryValues[e]);
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

  /* Initialize the free list. */
  if (listmatrix->numNonzeros < listmatrix->memNonzeros)
  {
    listmatrix->firstFreeNonzero = &listmatrix->nonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < listmatrix->memNonzeros - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i + 1];
    listmatrix->nonzeros[listmatrix->memNonzeros-1].right = NULL;
  }
  else
    listmatrix->firstFreeNonzero = NULL;

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8PrintDense(CMR* cmr, ListMat8* listmatrix, FILE* stream)
{
  assert(cmr);
  assert(listmatrix);
  assert(stream);

  int8_t* dense = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &dense, listmatrix->numColumns) );
  for (size_t column = 0; column < listmatrix->numColumns; ++column)
    dense[column] = 0;

  for (size_t row = 0; row < listmatrix->numRows; ++row)
  {
    for (ListMat8Nonzero* nz = listmatrix->rowElements[row].head.right;
      nz != &listmatrix->rowElements[row].head; nz = nz->right)
    {
      assert(nz->row == row);
      dense[nz->column] = nz->value;
    }
    for (size_t column = 0; column < listmatrix->numColumns; ++column)
    {
      fprintf(stream, " %3d", dense[column]);
      dense[column] = 0;
    }
    fprintf(stream, "\n");
  }
  fflush(stream);

  CMR_CALL( CMRfreeStackArray(cmr, &dense) );

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64PrintDense(CMR* cmr, ListMat64* listmatrix, FILE* stream)
{
  assert(cmr);
  assert(listmatrix);
  assert(stream);

  int64_t* dense = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &dense, listmatrix->numColumns) );
  for (size_t column = 0; column < listmatrix->numColumns; ++column)
    dense[column] = 0;

  for (size_t row = 0; row < listmatrix->numRows; ++row)
  {
    for (ListMat64Nonzero* nz = listmatrix->rowElements[row].head.right;
      nz != &listmatrix->rowElements[row].head; nz = nz->right)
    {
      assert(nz->row == row);
      dense[nz->column] = nz->value;
    }
    for (size_t column = 0; column < listmatrix->numColumns; ++column)
    {
      fprintf(stream, " %8ld", dense[column]);
      dense[column] = 0;
    }
    fprintf(stream, " (%ld nonzeros)\n", listmatrix->rowElements[row].numNonzeros);
  }
  fflush(stream);

  CMR_CALL( CMRfreeStackArray(cmr, &dense) );

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPPrintDense(CMR* cmr, ListMatGMP* listmatrix, FILE* stream)
{
  assert(cmr);
  assert(listmatrix);
  assert(stream);

  mpz_t* dense = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &dense, listmatrix->numColumns) );
  for (size_t column = 0; column < listmatrix->numColumns; ++column)
    mpz_init(dense[column]);

  for (size_t row = 0; row < listmatrix->numRows; ++row)
  {
    for (ListMatGMPNonzero* nz = listmatrix->rowElements[row].head.right;
      nz != &listmatrix->rowElements[row].head; nz = nz->right)
    {
      assert(nz->row == row);
      mpz_set(dense[nz->column], nz->value);
    }
    for (size_t column = 0; column < listmatrix->numColumns; ++column)
    {
      char buffer[1024];
      mpz_get_str(buffer, 10, dense[column]);
      fprintf(stream, "%6s", buffer);
      mpz_set_si(dense[column], 0);
    }
    fprintf(stream, "\n");
  }
  fflush(stream);

  for (size_t column = 0; column < listmatrix->numColumns; ++column)
    mpz_clear(dense[column]);
  CMR_CALL( CMRfreeStackArray(cmr, &dense) );

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8Insert(CMR* cmr, ListMat8* listmatrix, size_t row, size_t column, int8_t value, long special,
  ptrdiff_t* pmemoryShift)
{
  assert(cmr);
  assert(listmatrix);
  assert(row < listmatrix->numRows);
  assert(row < listmatrix->numColumns);

  if (!listmatrix->firstFreeNonzero)
  {
    assert(listmatrix->numNonzeros == listmatrix->memNonzeros);
    size_t newSize = 2 * listmatrix->memNonzeros;
    if (newSize < 256)
      newSize = 256;
        ListMat8Nonzero* newNonzeros = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &newNonzeros, newSize) );
    ptrdiff_t memoryShift = newNonzeros - listmatrix->nonzeros;

    /* Shift all pointers to other nonzeros by difference between new and old nonzero memory. */
    for (size_t i = 0; i < listmatrix->numNonzeros; ++i)
    {
      newNonzeros[i] = listmatrix->nonzeros[i];
      newNonzeros[i].left += memoryShift;
      newNonzeros[i].right += memoryShift;
      newNonzeros[i].above += memoryShift;
      newNonzeros[i].below += memoryShift;
    }

    /* Also do this for the heads, but repair it for their predecessors / successors. */
    for (size_t row = 0; row < listmatrix->numRows; ++row)
    {
      ListMat8Element* element = &listmatrix->rowElements[row];
      if (element->numNonzeros)
      {
        element->head.right += memoryShift;
        element->head.left += memoryShift;
        element->head.right->left -= memoryShift;
        element->head.left->right -= memoryShift;
      }
    }
    for (size_t column = 0; column < listmatrix->numColumns; ++column)
    {
      ListMat8Element* element = &listmatrix->columnElements[column];
      if (element->numNonzeros)
      {
        element->head.below += memoryShift;
        element->head.above += memoryShift;
        element->head.below->above -= memoryShift;
        element->head.above->below -= memoryShift;
      }
    }

    /* Create free list. */
    listmatrix->firstFreeNonzero = &newNonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < newSize - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i+1];
    listmatrix->nonzeros[newSize-1].right = NULL;

    /* Move old nonzero array. */
    listmatrix->memNonzeros = newSize;
    CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->nonzeros) );
    listmatrix->nonzeros = newNonzeros;

    if (pmemoryShift)
      *pmemoryShift = memoryShift;
  }
  else if (pmemoryShift)
    *pmemoryShift = 0;

  /* Actually insert the element. */
  ListMat8Nonzero* nz = listmatrix->firstFreeNonzero;
  assert(nz);
  listmatrix->firstFreeNonzero = nz->right;
  nz->row = row;
  nz->column = column;
  nz->value = value;
  nz->special = special;

  ListMat8Nonzero* head = &listmatrix->rowElements[row].head;
  nz->left = head;
  nz->right = head->right;
  head->right->left = nz;
  head->right = nz;

  head = &listmatrix->columnElements[column].head;
  nz->above = head;
  nz->below = head->below;
  head->below->above = nz;
  head->below = nz;

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64Insert(CMR* cmr, ListMat64* listmatrix, size_t row, size_t column, int64_t value,
  long special, ptrdiff_t* pmemoryShift)
{
  assert(cmr);
  assert(listmatrix);
  assert(row < listmatrix->numRows);
  assert(row < listmatrix->numColumns);

  CMRdbgMsg(10, "CMRlistmat64Insert for %ld of %ld possible nonzeros.\n", listmatrix->numNonzeros,
    listmatrix->memNonzeros);

  if (!listmatrix->firstFreeNonzero)
  {
    assert(listmatrix->numNonzeros == listmatrix->memNonzeros);
    size_t newSize = 2 * listmatrix->memNonzeros;
    if (newSize < 256)
      newSize = 256;
    ListMat64Nonzero* newNonzeros = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &newNonzeros, newSize) );
    ptrdiff_t memoryShift = newNonzeros - listmatrix->nonzeros;

    /* Shift all pointers to other nonzeros by difference between new and old nonzero memory. */
    for (size_t i = 0; i < listmatrix->numNonzeros; ++i)
    {
      newNonzeros[i] = listmatrix->nonzeros[i];
      newNonzeros[i].left += memoryShift;
      newNonzeros[i].right += memoryShift;
      newNonzeros[i].above += memoryShift;
      newNonzeros[i].below += memoryShift;
    }

    /* Also do this for the heads, but repair it for their predecessors / successors. */
    for (size_t row = 0; row < listmatrix->numRows; ++row)
    {
      ListMat64Element* element = &listmatrix->rowElements[row];
      if (element->numNonzeros)
      {
        element->head.right += memoryShift;
        element->head.left += memoryShift;
        element->head.right->left -= memoryShift;
        element->head.left->right -= memoryShift;
      }
    }
    for (size_t column = 0; column < listmatrix->numColumns; ++column)
    {
      ListMat64Element* element = &listmatrix->columnElements[column];
      if (element->numNonzeros)
      {
        element->head.below += memoryShift;
        element->head.above += memoryShift;
        element->head.below->above -= memoryShift;
        element->head.above->below -= memoryShift;
      }
    }

    /* Create free list. */
    listmatrix->firstFreeNonzero = &newNonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < newSize - 1; ++i)
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i+1];
    listmatrix->nonzeros[newSize-1].right = NULL;

    /* Move old nonzero array. */
    listmatrix->memNonzeros = newSize;
    CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->nonzeros) );
    listmatrix->nonzeros = newNonzeros;

    if (pmemoryShift)
      *pmemoryShift = memoryShift;
  }
  else if (pmemoryShift)
    *pmemoryShift = 0;

  /* Actually insert the element. */
  ListMat64Nonzero* nz = listmatrix->firstFreeNonzero;
  assert(nz);
  listmatrix->firstFreeNonzero = nz->right;
  nz->row = row;
  nz->column = column;
  nz->value = value;
  nz->special = special;

  ListMat64Nonzero* head = &listmatrix->rowElements[row].head;
  CMRdbgMsg(10, "Inserting %ld in column %ld horizontally between head and %ld in column %ld.\n", nz->value, nz->column,
    head->right->value, head->right->column);
  nz->left = head;
  nz->right = head->right;
  head->right->left = nz;
  head->right = nz;
  listmatrix->rowElements[row].numNonzeros++;

  head = &listmatrix->columnElements[column].head;
  CMRdbgMsg(10, "Inserting %ld in row %ld vertically between head and %ld in row %ld.\n", nz->value, nz->row,
    head->below->value, head->below->column);
  nz->above = head;
  nz->below = head->below;
  head->below->above = nz;
  head->below = nz;
  listmatrix->columnElements[column].numNonzeros++;

  listmatrix->numNonzeros++;

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPInsert(CMR* cmr, ListMatGMP* listmatrix, size_t row, size_t column, mpz_srcptr value,
  long special, ptrdiff_t* pmemoryShift)
{
  assert(cmr);
  assert(listmatrix);
  assert(row < listmatrix->numRows);
  assert(row < listmatrix->numColumns);

  CMRdbgMsg(10, "CMRlistmatGMPInsert for %ld of %ld possible nonzeros.\n", listmatrix->numNonzeros,
    listmatrix->memNonzeros);

  if (!listmatrix->firstFreeNonzero)
  {
    assert(listmatrix->numNonzeros == listmatrix->memNonzeros);
    size_t newSize = 2 * listmatrix->memNonzeros;
    if (newSize < 256)
      newSize = 256;
    ListMatGMPNonzero* newNonzeros = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &newNonzeros, newSize) );
    ptrdiff_t memoryShift = newNonzeros - listmatrix->nonzeros;

    /* Shift all pointers to other nonzeros by difference between new and old nonzero memory. */
    for (size_t i = 0; i < listmatrix->numNonzeros; ++i)
    {
      newNonzeros[i] = listmatrix->nonzeros[i];
      newNonzeros[i].left += memoryShift;
      newNonzeros[i].right += memoryShift;
      newNonzeros[i].above += memoryShift;
      newNonzeros[i].below += memoryShift;
      mpz_init(newNonzeros[i].value);
      mpz_swap(newNonzeros[i].value, listmatrix->nonzeros[i].value);
      mpz_clear(listmatrix->nonzeros[i].value);
    }

    /* Also do this for the heads, but repair it for their predecessors / successors. */
    for (size_t row = 0; row < listmatrix->numRows; ++row)
    {
      ListMatGMPElement* element = &listmatrix->rowElements[row];
      if (element->numNonzeros)
      {
        element->head.right += memoryShift;
        element->head.left += memoryShift;
        element->head.right->left -= memoryShift;
        element->head.left->right -= memoryShift;
      }
    }
    for (size_t column = 0; column < listmatrix->numColumns; ++column)
    {
      ListMatGMPElement* element = &listmatrix->columnElements[column];
      if (element->numNonzeros)
      {
        element->head.below += memoryShift;
        element->head.above += memoryShift;
        element->head.below->above -= memoryShift;
        element->head.above->below -= memoryShift;
      }
    }

    /* Create free list. */
    listmatrix->firstFreeNonzero = &newNonzeros[listmatrix->numNonzeros];
    for (size_t i = listmatrix->numNonzeros; i < newSize - 1; ++i)
    {
      listmatrix->nonzeros[i].right = &listmatrix->nonzeros[i+1];
      mpz_init(listmatrix->nonzeros[i].value);
    }
    listmatrix->nonzeros[newSize-1].right = NULL;
    mpz_init(listmatrix->nonzeros[newSize-1].value);

    /* Move old nonzero array. */
    listmatrix->memNonzeros = newSize;
    CMR_CALL( CMRfreeBlockArray(cmr, &listmatrix->nonzeros) );
    listmatrix->nonzeros = newNonzeros;

    if (pmemoryShift)
      *pmemoryShift = memoryShift;
  }
  else if (pmemoryShift)
    *pmemoryShift = 0;

  /* Actually insert the element. */
  ListMatGMPNonzero* nz = listmatrix->firstFreeNonzero;
  assert(nz);
  listmatrix->firstFreeNonzero = nz->right;
  nz->row = row;
  nz->column = column;
  mpz_set(nz->value, value);
  nz->special = special;

  ListMatGMPNonzero* head = &listmatrix->rowElements[row].head;
  CMRdbgMsg(10, "Inserting %ld in column %ld horizontally between head and %ld in column %ld.\n", nz->value, nz->column,
    head->right->value, head->right->column);
  nz->left = head;
  nz->right = head->right;
  head->right->left = nz;
  head->right = nz;
  listmatrix->rowElements[row].numNonzeros++;

  head = &listmatrix->columnElements[column].head;
  CMRdbgMsg(10, "Inserting %ld in row %ld vertically between head and %ld in row %ld.\n", nz->value, nz->row,
    head->below->value, head->below->column);
  nz->above = head;
  nz->below = head->below;
  head->below->above = nz;
  head->below = nz;
  listmatrix->columnElements[column].numNonzeros++;

  listmatrix->numNonzeros++;

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

CMR_ERROR CMRlistmat8Delete(CMR* cmr, ListMat8* listmatrix, ListMat8Nonzero* nz)
{
  assert(cmr);
  assert(listmatrix);
  assert(nz);
  assert(nz->row < listmatrix->numRows);
  assert(nz->column < listmatrix->numColumns);
  assert(nz != &listmatrix->rowElements[nz->row].head);
  assert(nz != &listmatrix->columnElements[nz->column].head);

  listmatrix->numNonzeros--;
  listmatrix->rowElements[nz->row].numNonzeros--;
  listmatrix->columnElements[nz->column].numNonzeros--;

  nz->left->right = nz->right;
  nz->right->left = nz->left;
  nz->above->below = nz->below;
  nz->below->above = nz->above;
  nz->right = listmatrix->firstFreeNonzero;
  listmatrix->firstFreeNonzero = nz;
  listmatrix->numNonzeros--;

  return CMR_OKAY;
}

CMR_ERROR CMRlistmat64Delete(CMR* cmr, ListMat64* listmatrix, ListMat64Nonzero* nz)
{
  assert(cmr);
  assert(listmatrix);
  assert(nz);
  assert(nz->row < listmatrix->numRows);
  assert(nz->column < listmatrix->numColumns);
  assert(nz != &listmatrix->rowElements[nz->row].head);
  assert(nz != &listmatrix->columnElements[nz->column].head);

  listmatrix->numNonzeros--;
  listmatrix->rowElements[nz->row].numNonzeros--;
  listmatrix->columnElements[nz->column].numNonzeros--;

  nz->left->right = nz->right;
  nz->right->left = nz->left;
  nz->above->below = nz->below;
  nz->below->above = nz->above;
  nz->right = listmatrix->firstFreeNonzero;
  listmatrix->firstFreeNonzero = nz;

  return CMR_OKAY;
}

#if defined(CMR_WITH_GMP)

CMR_ERROR CMRlistmatGMPDelete(CMR* cmr, ListMatGMP* listmatrix, ListMatGMPNonzero* nz)
{
  assert(cmr);
  assert(listmatrix);
  assert(nz);
  assert(nz->row < listmatrix->numRows);
  assert(nz->column < listmatrix->numColumns);
  assert(nz != &listmatrix->rowElements[nz->row].head);
  assert(nz != &listmatrix->columnElements[nz->column].head);

  listmatrix->numNonzeros--;
  listmatrix->rowElements[nz->row].numNonzeros--;
  listmatrix->columnElements[nz->column].numNonzeros--;

  nz->left->right = nz->right;
  nz->right->left = nz->left;
  nz->above->below = nz->below;
  nz->below->above = nz->above;
  nz->right = listmatrix->firstFreeNonzero;
  listmatrix->firstFreeNonzero = nz;

  return CMR_OKAY;
}

#endif /* CMR_WITH_GMP */

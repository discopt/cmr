// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/separation.h>

#include "env_internal.h"

#include <stdint.h>
#include <assert.h>

CMR_ERROR CMRsepaCreate(CMR* cmr, size_t numRows, size_t numColumns, CMR_SEPA** psepa)
{
  assert(cmr);
  assert(psepa);

  size_t numElements = numRows + numColumns;

  CMR_CALL( CMRallocBlock(cmr, psepa) );
  CMR_SEPA* sepa = *psepa;
  sepa->indicatorMemory = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &sepa->indicatorMemory, numElements) );
  sepa->elementMemory = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &sepa->elementMemory, numElements) );
  sepa->numRows[0] = 0;
  sepa->numRows[1] = numRows;
  sepa->numColumns[0] = 0;
  sepa->numColumns[1] = numColumns;
  sepa->rowsToPart = &sepa->indicatorMemory[0];
  sepa->columnsToPart = &sepa->indicatorMemory[numRows];

  return CMR_OKAY;
}

CMR_ERROR CMRsepaInitialize(CMR* cmr, CMR_SEPA* sepa, size_t firstExtraRow0, size_t firstExtraColumn1,
  size_t firstExtraRow1, size_t firstExtraColumn0, size_t secondExtraRow0, size_t secondExtraColumn1,
  size_t secondExtraRow1, size_t secondExtraColumn0)
{
  assert(cmr);
  assert(sepa);
  assert((firstExtraRow0 < SIZE_MAX && firstExtraColumn1 < SIZE_MAX)
    || (firstExtraRow0 == SIZE_MAX && firstExtraColumn1 == SIZE_MAX));
  assert((firstExtraRow1 < SIZE_MAX && firstExtraColumn0 < SIZE_MAX)
    || (firstExtraRow1 == SIZE_MAX && firstExtraColumn0 == SIZE_MAX));
  assert((secondExtraRow0 < SIZE_MAX && secondExtraColumn1 < SIZE_MAX)
    || (secondExtraRow0 == SIZE_MAX && secondExtraColumn1 == SIZE_MAX));
  assert((secondExtraRow1 < SIZE_MAX && secondExtraColumn0 < SIZE_MAX)
    || (secondExtraRow1 == SIZE_MAX && secondExtraColumn0 == SIZE_MAX));

  size_t numRows = sepa->numRows[0] + sepa->numRows[1];
  size_t numColumns = sepa->numColumns[0] + sepa->numColumns[1];

  /* Count sizes of parts. */
  sepa->numRows[0] = 0;
  sepa->numRows[1] = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    unsigned char part = sepa->rowsToPart[row];
    if (part >= 2)
      continue;
    sepa->numRows[part]++;
  }
  sepa->numColumns[0] = 0;
  sepa->numColumns[1] = 0;
  for (size_t column = 0; column < numColumns; ++column)
  {
    unsigned char part = sepa->columnsToPart[column];
    if (part >= 2)
      continue;
    sepa->numColumns[part]++;
  }

  /* Set array pointers. */
  sepa->rows[0] = &sepa->elementMemory[0];
  sepa->rows[1] = &sepa->elementMemory[sepa->numRows[0]];
  sepa->columns[0] = &sepa->elementMemory[numRows];
  sepa->columns[1] = &sepa->elementMemory[numRows + sepa->numColumns[0]];

  /* Fill arrays. */
  sepa->numRows[0] = 0;
  sepa->numRows[1] = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    unsigned char part = sepa->rowsToPart[row];
    if (part >= 2)
      continue;
    sepa->rows[part][sepa->numRows[part]] = row;
    sepa->numRows[part]++;
  }
  sepa->numColumns[0] = 0;
  sepa->numColumns[1] = 0;
  for (size_t column = 0; column < numColumns; ++column)
  {
    unsigned char part = sepa->columnsToPart[column];
    if (part >= 2)
      continue;
    sepa->columns[part][sepa->numColumns[part]] = column;
    sepa->numColumns[part]++;
  }

  sepa->extraRows0[0] = firstExtraRow0;
  sepa->extraColumns1[0] = firstExtraColumn1;
  sepa->extraRows1[0] = firstExtraRow1;
  sepa->extraColumns0[0] = firstExtraColumn0;
  sepa->extraRows0[1] = secondExtraRow0;
  sepa->extraColumns1[1] = secondExtraColumn1;
  sepa->extraRows1[1] = secondExtraRow1;
  sepa->extraColumns0[1] = secondExtraColumn0;

  return CMR_OKAY;
}

CMR_ERROR CMRsepaFree(CMR* cmr, CMR_SEPA** psepa)
{
  assert(cmr);
  assert(psepa);

  if (*psepa)
  {
    CMR_SEPA* sepa = *psepa;
    CMR_CALL( CMRfreeBlockArray(cmr, &sepa->elementMemory) );
    CMR_CALL( CMRfreeBlockArray(cmr, &sepa->indicatorMemory) );
    CMR_CALL( CMRfreeBlock(cmr, psepa) );
  }

  return CMR_OKAY;
}

static
CMR_ERROR checkTernary(
  CMR* cmr,
  size_t numBlockRows,
  size_t* blockRows,
  size_t numBlockColumns,
  size_t* blockColumns,
  unsigned char rank,
  CMR_CHRMAT* matrix,
  CMR_SUBMAT* submatrix,
  bool* pisTernary,
  CMR_SUBMAT** psubmatrix
)
{
  assert(cmr);
  assert(blockRows);
  assert(blockColumns);
  assert(rank <= 3);
  assert(matrix);
  assert(pisTernary);

  *pisTernary = true;
  if (rank == 0)
    return CMR_OKAY;

  if (rank == 1)
  {
    CMRdbgMsg(2, "Checking a binary rank-1 %ldx%ld submatrix for ternary rank 1.\n", numBlockRows, numBlockColumns);

    size_t* entries = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &entries, numBlockRows) );
    for (size_t i = 0; i < numBlockRows; ++i)
    {
      size_t blockRow = submatrix ? submatrix->rows[blockRows[i]] : blockRows[i];
      entries[i] = matrix->rowSlice[blockRow];
    }

    size_t firstRow = SIZE_MAX;
    size_t firstColumn = SIZE_MAX;
    char* firstColumnValues = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &firstColumnValues, numBlockRows) );
    for (size_t i = 0; i < numBlockRows; ++i)
      firstColumnValues[i] = 0;

    for (size_t j = 0; j < numBlockColumns && *pisTernary; ++j)
    {
      size_t column = submatrix ? submatrix->columns[blockColumns[j]] : blockColumns[j];
      bool negated;
      for (size_t i = 0; i < numBlockRows && *pisTernary; ++i)
      {
        size_t row = blockRows[i];
        while (matrix->entryColumns[entries[i]] < column)
          ++entries[i];
        if (matrix->entryColumns[entries[i]] == column)
        {
          if (firstColumn == SIZE_MAX)
          {
            firstRow = row;
            firstColumn = column;
          }

          if (column == firstColumn)
          {
            firstColumnValues[i] = matrix->entryValues[entries[i]];
          }
          else if (row == firstRow)
            negated = firstColumnValues[i] != matrix->entryValues[entries[i]];
          else
          {
            bool neg = firstColumnValues[i] != matrix->entryValues[entries[i]];
            if (negated != neg)
            {
              CMR_CALL( CMRsubmatCreate(cmr, 2, 2, psubmatrix) );
              CMR_SUBMAT* submatrix = *psubmatrix;
              submatrix->rows[0] = firstRow;
              submatrix->rows[1] = row;
              submatrix->columns[0] = firstColumn;
              submatrix->columns[1] = column;
              *pisTernary = false;
            }
          }
        }
      }
    }
    
    CMR_CALL( CMRfreeStackArray(cmr, &firstColumnValues) );
    CMR_CALL( CMRfreeStackArray(cmr, &entries) );
  }
  else
  {
    assert("Checking separation for being ternary is not implemented for rank 2." == 0);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRsepaCheckTernary(CMR* cmr, CMR_SEPA* sepa, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix, bool* pisTernary,
  CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(sepa);
  assert(matrix);
  assert(pisTernary);
  assert(!psubmatrix || !*psubmatrix);

  CMR_CALL( checkTernary(cmr, sepa->numRows[0], sepa->rows[0], sepa->numColumns[1], sepa->columns[1],
    CMRsepaRankTopRight(sepa), matrix, submatrix, pisTernary, psubmatrix) );
  if (*pisTernary)
  {
    CMR_CALL( checkTernary(cmr, sepa->numRows[1], sepa->rows[1], sepa->numColumns[0], sepa->columns[0],
      CMRsepaRankBottomLeft(sepa), matrix, submatrix, pisTernary, psubmatrix) );
  }

  return CMR_OKAY;
}

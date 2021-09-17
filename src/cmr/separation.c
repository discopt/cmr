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

static
CMR_ERROR initialize(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_SEPA* separation  /**< Already created separation. */
)
{
  assert(cmr);
  assert(separation);

  size_t numRows = separation->numRows[0] + separation->numRows[1];
  size_t numColumns = separation->numColumns[0] + separation->numColumns[1];

  /* Count sizes of parts. */
  separation->numRows[0] = 0;
  separation->numRows[1] = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    unsigned char part = separation->rowsToPart[row];
    if (part >= 2)
      continue;
    separation->numRows[part]++;
  }
  separation->numColumns[0] = 0;
  separation->numColumns[1] = 0;
  for (size_t column = 0; column < numColumns; ++column)
  {
    unsigned char part = separation->columnsToPart[column];
    if (part >= 2)
      continue;
    separation->numColumns[part]++;
  }

  /* Set array pointers. */
  separation->rows[0] = &separation->elementMemory[0];
  separation->rows[1] = &separation->elementMemory[separation->numRows[0]];
  separation->columns[0] = &separation->elementMemory[numRows];
  separation->columns[1] = &separation->elementMemory[numRows + separation->numColumns[0]];

  /* Fill arrays. */
  separation->numRows[0] = 0;
  separation->numRows[1] = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    unsigned char part = separation->rowsToPart[row];
    if (part >= 2)
      continue;
    separation->rows[part][separation->numRows[part]] = row;
    separation->numRows[part]++;
  }
  separation->numColumns[0] = 0;
  separation->numColumns[1] = 0;
  for (size_t column = 0; column < numColumns; ++column)
  {
    unsigned char part = separation->columnsToPart[column];
    if (part >= 2)
      continue;
    separation->columns[part][separation->numColumns[part]] = column;
    separation->numColumns[part]++;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRsepaInitialize(CMR* cmr, CMR_SEPA* separation, size_t firstExtraRow0, size_t firstExtraColumn1,
  size_t firstExtraRow1, size_t firstExtraColumn0, size_t secondExtraRow0, size_t secondExtraColumn1,
  size_t secondExtraRow1, size_t secondExtraColumn0)
{
  assert(cmr);
  assert(separation);
  assert((firstExtraRow0 < SIZE_MAX && firstExtraColumn1 < SIZE_MAX)
    || (firstExtraRow0 == SIZE_MAX && firstExtraColumn1 == SIZE_MAX));
  assert((firstExtraRow1 < SIZE_MAX && firstExtraColumn0 < SIZE_MAX)
    || (firstExtraRow1 == SIZE_MAX && firstExtraColumn0 == SIZE_MAX));
  assert((secondExtraRow0 < SIZE_MAX && secondExtraColumn1 < SIZE_MAX)
    || (secondExtraRow0 == SIZE_MAX && secondExtraColumn1 == SIZE_MAX));
  assert((secondExtraRow1 < SIZE_MAX && secondExtraColumn0 < SIZE_MAX)
    || (secondExtraRow1 == SIZE_MAX && secondExtraColumn0 == SIZE_MAX));

  CMR_CALL( initialize(cmr, separation) );

  separation->extraRows0[0] = firstExtraRow0;
  separation->extraColumns1[0] = firstExtraColumn1;
  separation->extraRows1[0] = firstExtraRow1;
  separation->extraColumns0[0] = firstExtraColumn0;
  separation->extraRows0[1] = secondExtraRow0;
  separation->extraColumns1[1] = secondExtraColumn1;
  separation->extraRows1[1] = secondExtraRow1;
  separation->extraColumns0[1] = secondExtraColumn0;

  return CMR_OKAY;
}

CMR_ERROR CMRsepaInitializeMatrix(CMR* cmr, CMR_SEPA* separation, CMR_CHRMAT* matrix, unsigned char totalRank)
{
  assert(cmr);
  assert(separation);
  assert(matrix);

  CMR_CALL( initialize(cmr, separation) );

  if (totalRank == 0)
    return CMR_OKAY;

  assert(totalRank <= 1 || "Not implemented" == 0);

  CMRdbgMsg(4, "Searching for rank-1 submatrix in matrix\n");
#if defined(CMR_DEBUG)
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  for (int i = 0; i < 2; ++i)
  {
    separation->extraRows0[i] = SIZE_MAX;
    separation->extraRows1[i] = SIZE_MAX;
    separation->extraColumns0[i] = SIZE_MAX;
    separation->extraColumns1[i] = SIZE_MAX;
  }

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    unsigned char rowPart = separation->rowsToPart[row];
    CMRdbgMsg(6, "Row r%ld is in part %d\n", row+1, rowPart);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      CMRdbgMsg(6, "Column c%ld is in part %d\n", column+1, separation->columnsToPart[column]);
      if (rowPart != separation->columnsToPart[column])
      {
        if (rowPart == 0)
        {
          separation->extraColumns0[0] = column;
          separation->extraRows1[0] = row;
        }
        else
        {
          separation->extraRows0[0] = row;
          separation->extraColumns1[0] = column;
        }
        return CMR_OKAY;
      }
    }
  }

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

CMR_ERROR CMRoneSum(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, CMR_CHRMAT ** presult)
{
  assert(cmr);
  assert(first);
  assert(second);
  assert(presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, first->numRows + second->numRows, first->numColumns + second->numColumns,
    first->numNonzeros + second->numNonzeros) );

  CMR_CHRMAT* result = *presult;

  for (size_t row = 0; row < first->numRows; ++row)
    result->rowSlice[row] = first->rowSlice[row];
  for (size_t row = 0; row <= second->numRows; ++row)
    result->rowSlice[first->numRows + row] = second->rowSlice[row] + first->numNonzeros;

  for (size_t e = 0; e < first->numNonzeros; ++e)
  {
    result->entryColumns[e] = first->entryColumns[e];
    result->entryValues[e] = first->entryValues[e];
  }

  for (size_t e = 0; e < second->numNonzeros; ++e)
  {
    result->entryColumns[first->numNonzeros + e] = first->numColumns + second->entryColumns[e];
    result->entryValues[first->numNonzeros + e] = second->entryValues[e];
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtwoSum(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, CMR_ELEMENT firstMarker, CMR_ELEMENT secondMarker, CMR_CHRMAT ** presult)
{
  assert(cmr);
  assert(first);
  assert(second);
  assert(presult);

  if ((CMRelementIsRow(firstMarker) && CMRelementIsRow(secondMarker))
    || (CMRelementIsColumn(firstMarker) && CMRelementIsColumn(secondMarker)))
  {
    return CMR_ERROR_INPUT;
  }

  bool bottomLeft = CMRelementIsRow(firstMarker);
  size_t firstRowMarker = SIZE_MAX;
  size_t secondRowMarker = SIZE_MAX;
  size_t firstColumnMarker = SIZE_MAX;
  size_t secondColumnMarker = SIZE_MAX;
  char* markerColumn = NULL;
  size_t markerColumnNumNonzeros = 0;
  size_t markerRowNumNonzeros = 0;
  if (bottomLeft)
  {
    firstRowMarker = CMRelementToRowIndex(firstMarker);
    markerRowNumNonzeros = first->rowSlice[firstRowMarker+1] - first->rowSlice[firstRowMarker];

    secondColumnMarker = CMRelementToColumnIndex(secondMarker);
    CMR_CALL( CMRallocStackArray(cmr, &markerColumn, second->numColumns) );
    for (size_t row = 0; row < second->numRows; ++row)
    {
      size_t entry;
      CMR_CALL( CMRchrmatFindEntry(cmr, second, row, secondColumnMarker, &entry) );
      if (entry == SIZE_MAX)
      {
        markerColumn[row] = 0;
        CMRdbgMsg(2, "markerColumn[%ld] = zero.\n", row);
      }
      else
      {
        assert(second->entryColumns[entry] == secondColumnMarker);
        markerColumn[row] = second->entryValues[entry];
        CMRdbgMsg(2, "markerColumn[%ld] = %d.\n", row, second->entryValues[entry]);
        ++markerColumnNumNonzeros;
      }
    }
  }
  else
  {
    firstColumnMarker = CMRelementToColumnIndex(firstMarker); 
    CMR_CALL( CMRallocStackArray(cmr, &markerColumn, first->numColumns) );
    for (size_t row = 0; row < first->numRows; ++row)
    {
      size_t entry;
      CMR_CALL( CMRchrmatFindEntry(cmr, first, row, firstColumnMarker, &entry) );
      if (entry == SIZE_MAX)
        markerColumn[row] = 0;
      else
      {
        assert(first->entryColumns[entry] == firstColumnMarker);
        markerColumn[row] = first->entryValues[entry];
        ++markerColumnNumNonzeros;
      }
    }

    secondRowMarker = CMRelementToRowIndex(secondMarker);
    markerRowNumNonzeros = second->rowSlice[secondRowMarker+1] - second->rowSlice[secondRowMarker];
  }

  /* Create the actual matrix. */
  CMR_CALL( CMRchrmatCreate(cmr, presult, first->numRows + second->numRows - 1,
    first->numColumns + second->numColumns - 1,
    first->numNonzeros + second->numNonzeros + (markerRowNumNonzeros - 1) * (markerColumnNumNonzeros - 1) - 1) );
  CMR_CHRMAT* result = *presult;
  size_t resultRow = 0;
  size_t resultNonzero = 0;
  for (size_t row = 0; row < first->numRows; ++row)
  {
    if (row != firstRowMarker)
    {
      result->rowSlice[resultRow] = resultNonzero;
      CMRdbgMsg(0, "Row %ld starts at nonzero #%ld.\n", resultRow, resultNonzero);

      /* First matrix is top-left. */
      for (size_t e = first->rowSlice[row]; e < first->rowSlice[row + 1]; ++e)
      {
        size_t column = first->entryColumns[e];
        if (column == firstColumnMarker)
          continue;
        result->entryValues[resultNonzero] = first->entryValues[e];
        result->entryColumns[resultNonzero] = column < firstColumnMarker ? column : (column - 1);
        CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
        ++resultNonzero;
      }

      /* Rank-1 matrix is top-right. */
      if (!bottomLeft && markerColumn[row])
      {
        for (size_t e = second->rowSlice[secondRowMarker]; e < second->rowSlice[secondRowMarker + 1]; ++e)
        {
          result->entryValues[resultNonzero] = markerColumn[row] * second->entryValues[e];
          result->entryColumns[resultNonzero] = first->numColumns - 1 + second->entryColumns[e];
          CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
          ++resultNonzero;
        }
      }
      ++resultRow;
    }
  }
  for (size_t row = 0; row < second->numRows; ++row)
  {
    if (row != secondRowMarker)
    {
      result->rowSlice[resultRow] = resultNonzero;
      CMRdbgMsg(0, "Row %ld starts at nonzero #%ld.\n", resultRow, resultNonzero);

      /* Rank-1 matrix is bottom-left. */
      if (bottomLeft && markerColumn[row])
      {
        for (size_t e = first->rowSlice[firstRowMarker]; e < first->rowSlice[firstRowMarker + 1]; ++e)
        {
          size_t column = first->entryColumns[e];
          result->entryValues[resultNonzero] = markerColumn[row] * first->entryValues[e];
          result->entryColumns[resultNonzero] = column;
          CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
          ++resultNonzero;
        }
      }
      
      /* Second matrix is bottom-right. */
      for (size_t e = second->rowSlice[row]; e < second->rowSlice[row + 1]; ++e)
      {
        size_t column = second->entryColumns[e];
        if (column == secondColumnMarker)
          continue;
        result->entryValues[resultNonzero] = second->entryValues[e];
        result->entryColumns[resultNonzero] = (bottomLeft ? first->numColumns : (first->numColumns - 1))
           + (column < secondColumnMarker ? column : (column - 1));
        CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
        ++resultNonzero;
      }

      ++resultRow;
    }
  }
  result->rowSlice[result->numRows] = resultNonzero;
  assert(resultNonzero == result->numNonzeros);
  
  CMRconsistencyAssert( CMRchrmatConsistency(result) );

  CMR_CALL( CMRfreeStackArray(cmr, &markerColumn) );

  return CMR_OKAY;
}

CMR_ERROR CMRthreeSum(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, CMR_ELEMENT firstMarker1, CMR_ELEMENT secondMarker1, CMR_ELEMENT firstMarker2, CMR_ELEMENT secondMarker2, CMR_CHRMAT ** presult)
{
  assert(cmr);
  assert(first);
  assert(second);
  assert(presult);

  if ((CMRelementIsRow(firstMarker1) && CMRelementIsRow(secondMarker1))
    || (CMRelementIsColumn(firstMarker1) && CMRelementIsColumn(secondMarker1))
    || (CMRelementIsRow(firstMarker2) && CMRelementIsRow(secondMarker2))
    || (CMRelementIsColumn(firstMarker2) && CMRelementIsColumn(secondMarker2)))
  {
    return CMR_ERROR_INPUT;
  }

  assert("3-sums are not implemented, yet." == 0);

  return CMR_OKAY;
}



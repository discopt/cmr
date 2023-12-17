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
  CMR_UNUSED(cmr);

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

  separation->extraRows[0][0] = firstExtraRow0;
  separation->extraColumns[1][0] = firstExtraColumn1;
  separation->extraRows[1][0] = firstExtraRow1;
  separation->extraColumns[0][0] = firstExtraColumn0;
  separation->extraRows[0][1] = secondExtraRow0;
  separation->extraColumns[1][1] = secondExtraColumn1;
  separation->extraRows[1][1] = secondExtraRow1;
  separation->extraColumns[0][1] = secondExtraColumn0;

  return CMR_OKAY;
}

/**
 * \brief Checks whether the submatrix with rows assigned to \p part and columns assigned to the other has at least
 *        rank 1.
 * 
 * If the rank is at least 1, the row and column of a corresponding nonzero are stored in \p separation.
 */

static
bool findRank1(
  CMR_CHRMAT* matrix,     /**< Matrix. */
  CMR_SEPA* separation,   /**< Separation. */
  short part              /**< Part to which the investigated rows belong. */
)
{
  assert(matrix);
  assert(part >=0 && part < 2);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (separation->rowsToPart[row] == part)
      continue;
    size_t entry = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    while (entry < beyond)
    {
      size_t column = matrix->entryColumns[entry];
      if (separation->columnsToPart[column] == part)
      {
        separation->extraRows[part][0] = row;
        separation->extraColumns[1-part][0] = column;
        return true;
      }
      ++entry;
    }
  }

  return false;
}

/**
 * \brief Checks whether the submatrix with rows assigned to \p part and columns assigned to the other has at least
 *        rank 2.
 * 
 * Assumes that \p separation already contains a row and a column of some nonzero.
 * If the rank is at least 2, the row and column of an nonzero with which the previous one induces a rank-2 submatrix
 * are stored in \p separation.
 */

static
bool findRank2(
  CMR_CHRMAT* matrix,     /**< Matrix. */
  CMR_SEPA* separation,   /**< Separation. */
  short part              /**< Part to which the investigated rows belong. */
)
{
  assert(matrix);
  assert(separation);
  assert(part >=0 && part < 2);

  for (size_t row = separation->extraRows[part][0] + 1; row < matrix->numRows; ++row)
  {
    if (separation->rowsToPart[row] == part)
      continue;

    size_t entry = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    size_t column = (entry < beyond) ? matrix->entryColumns[entry] : SIZE_MAX;
    size_t entryRep = matrix->rowSlice[separation->extraRows[part][0]];
    size_t beyondRep = matrix->rowSlice[separation->extraRows[part][0] + 1];
    size_t columnRep = (entryRep < beyondRep) ? matrix->entryColumns[entryRep] : SIZE_MAX;
    bool isZero = true;
    bool equalRep = true;
    while (column != SIZE_MAX || columnRep != SIZE_MAX)
    {
      if (column < columnRep)
      {
        if (separation->columnsToPart[column] == part)
        {
          /* New row has a 1-entry but representative does not. */
          separation->extraRows[part][1] = row;
          separation->extraColumns[1-part][1] = column;
          return true;
        }
        else
        {
          /* We're outside of the submatrix. */
          ++entry;
          column = (entry < beyond) ? matrix->entryColumns[entry] : SIZE_MAX;
        }
      }
      else if (columnRep < column)
      {
        if (separation->columnsToPart[columnRep] == part)
        {
          if (!isZero)
          {
            /* We had a common nonzero before, so the 0 here yields rank 2. */
            separation->extraRows[part][1] = row;
            separation->extraColumns[1-part][1] = columnRep;
            return true;
          }
          else
            equalRep = false;
        }

        /* We're outside of the submatrix or still a zero vector. */
        ++entryRep;
        columnRep = (entryRep < beyondRep) ? matrix->entryColumns[entryRep] : SIZE_MAX;
      }
      else
      {
        if (separation->columnsToPart[column] == part)
        {
          if (!equalRep)
          {
            /* We had a 1 in rep with a 0 here before, so the two 1s here yield rank 2. */
            separation->extraRows[part][1] = row;
            separation->extraColumns[1-part][1] = columnRep;
            return true;
          }
          else
            isZero = false;
        }

        /* We're outside of the submatrix. */
        ++entry;
        column = (entry < beyond) ? matrix->entryColumns[entry] : SIZE_MAX;
        ++entryRep;
        columnRep = (entryRep < beyondRep) ? matrix->entryColumns[entryRep] : SIZE_MAX;
      }
    }
  }

  return false;
}

CMR_ERROR CMRsepaInitializeMatrix(CMR* cmr, CMR_SEPA* separation, CMR_CHRMAT* matrix, unsigned char totalRank)
{
  assert(cmr);
  assert(separation);
  assert(matrix);

  CMR_CALL( initialize(cmr, separation) );
  separation->extraRows[0][0] = SIZE_MAX;
  separation->extraRows[0][1] = SIZE_MAX;
  separation->extraRows[1][0] = SIZE_MAX;
  separation->extraRows[1][1] = SIZE_MAX;
  separation->extraColumns[0][0] = SIZE_MAX;
  separation->extraColumns[0][1] = SIZE_MAX;
  separation->extraColumns[1][0] = SIZE_MAX;
  separation->extraColumns[1][1] = SIZE_MAX;
  unsigned char ranks[2] = { 0, 0 };

#ifdef NDEBUG
  if (totalRank == 0)
    return CMR_OKAY;
#endif /* In debug mode we verify the rank. */

  /* Check submatrices for rank (at least) 1. */
  for (short part = 0; part < 2; ++part)
  {
    if (findRank1(matrix, separation, part))
      ++ranks[part];
#ifdef NDEBUG
    if (ranks[0] + ranks[1] == totalRank)
      return CMR_OKAY;
#endif
  }

  /* Check submatrices for rank (at least) 2. */
  for (short part = 0; part < 2; ++part)
  {
    if (ranks[part] == 0)
      continue;

    if (findRank2(matrix, separation, part))
      ++ranks[part];
#ifdef NDEBUG
    if (ranks[0] + ranks[1] == totalRank)
      return CMR_OKAY;
#endif
  }  

  if (ranks[0] + ranks[1] == totalRank)
    return CMR_OKAY;
  else
  {
    printf("Submatrix ranks are %d and %d, but expected total rank is %d.\n", ranks[0], ranks[1], totalRank);
    return CMR_ERROR_INVALID;
  }
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

CMR_ERROR CMRtwoSum(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, CMR_ELEMENT firstMarker, CMR_ELEMENT secondMarker,
  int8_t characteristic, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(first);
  assert(second);
  assert(presult);

  CMRdbgMsg(0, "CMRtwoSum of a %zux%zu and a %zux%zu matrix.\n", first->numRows, first->numColumns,
    second->numRows, second->numColumns);

  if ((CMRelementIsRow(firstMarker) && CMRelementIsRow(secondMarker))
    || (CMRelementIsColumn(firstMarker) && CMRelementIsColumn(secondMarker)))
  {
    return CMR_ERROR_INPUT;
  }

  bool bottomLeft = CMRelementIsRow(firstMarker);
  size_t firstRowMarker = SIZE_MAX; /* row index of A where to find a (or SIZE_MAX). */
  size_t secondRowMarker = SIZE_MAX; /* row index of B where to find b (or SIZE_MAX). */
  size_t firstColumnMarker = SIZE_MAX; /* column index of A where to find a (or SIZE_MAX). */
  size_t secondColumnMarker = SIZE_MAX; /* column index of B where to find b (or SIZE_MAX). */
  char* markerColumn = NULL; /* Nonzero entries of the column vector among a,b. */
  size_t markerColumnNumNonzeros = 0; /* Number of nonzeros in markerColumn. */
  size_t markerRowNumNonzeros = 0; /* Number of nonzeros of a,b that is not markerColumn. */
  CMR_ERROR error = CMR_OKAY;

  /* Step 1: fill marker column. */
  if (bottomLeft)
  {
    /* rank 1 in bottom-right matrix. */
    firstRowMarker = CMRelementToRowIndex(firstMarker);
    markerRowNumNonzeros = first->rowSlice[firstRowMarker+1] - first->rowSlice[firstRowMarker];

    secondColumnMarker = CMRelementToColumnIndex(secondMarker);
    CMR_CALL( CMRallocStackArray(cmr, &markerColumn, second->numColumns) );
    for (size_t row = 0; row < second->numRows; ++row)
    {
      size_t entry;
      CMR_CALL( CMRchrmatFindEntry(second, row, secondColumnMarker, &entry) );
      if (entry == SIZE_MAX)
      {
        markerColumn[row] = 0;
        CMRdbgMsg(2, "markerColumn[%ld] is 0.\n", row);
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
    /* rank 1 in top right */

    firstColumnMarker = CMRelementToColumnIndex(firstMarker); 
    CMR_CALL( CMRallocStackArray(cmr, &markerColumn, first->numColumns) );
    for (size_t row = 0; row < first->numRows; ++row)
    {
      size_t entry;
      CMR_CALL( CMRchrmatFindEntry(first, row, firstColumnMarker, &entry) );
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

  /* Step 2: create the actual matrix. */
  CMRdbgMsg(0, "#nonzeros is %zu + %zu + %zu * %zu - 1\n", first->numNonzeros, second->numNonzeros,
    markerRowNumNonzeros - 1, markerColumnNumNonzeros - 1);
  CMR_CALL( CMRchrmatCreate(cmr, presult, first->numRows + second->numRows - 1,
    first->numColumns + second->numColumns - 1,
    first->numNonzeros + second->numNonzeros + (markerRowNumNonzeros - 1) * (markerColumnNumNonzeros - 1) - 1) );
  CMR_CHRMAT* result = *presult;
  size_t resultRow = 0;
  size_t resultNonzero = 0;

  /* Step 2a: top part. */
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
          int resultValue = (int)(markerColumn[row]) * (int)second->entryValues[e];

          /* Take care of modulo operations. */
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }

          if (resultValue > INT8_MAX || resultValue < INT8_MIN)
          {
            error = CMR_ERROR_OVERFLOW;
            goto cleanup;
          }

          result->entryValues[resultNonzero] = resultValue;
          result->entryColumns[resultNonzero] = first->numColumns - 1 + second->entryColumns[e];
          CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
          ++resultNonzero;
        }
      }
      ++resultRow;
    }
  }

  /* Step 2b: bottom part. */
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

          /* Take care of modulo operations. */
          int resultValue = (int)markerColumn[row] * (int)first->entryValues[e];
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }

          if (resultValue > INT8_MAX || resultValue < INT8_MIN)
          {
            error = CMR_ERROR_OVERFLOW;
            goto cleanup;
          }

          result->entryValues[resultNonzero] = resultValue;
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

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &markerColumn) );

  return error;
}

CMR_ERROR CMRthreeSum(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, CMR_ELEMENT firstMarker1,
  CMR_ELEMENT secondMarker1, CMR_ELEMENT firstMarker2, CMR_ELEMENT secondMarker2, int8_t characteristic,
  CMR_CHRMAT** presult)
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

  if (CMRelementIsRow(firstMarker1) && CMRelementIsRow(firstMarker2)
    && (CMRelementToRowIndex(firstMarker1) == CMRelementToRowIndex(firstMarker2)))
  {
    return CMR_ERROR_INPUT;
  }
  if (CMRelementIsColumn(firstMarker1) && CMRelementIsColumn(firstMarker2)
    && (CMRelementToColumnIndex(firstMarker1) == CMRelementToColumnIndex(firstMarker2)))
  {
    return CMR_ERROR_INPUT;
  }
  if (CMRelementIsRow(secondMarker1) && CMRelementIsRow(secondMarker2)
    && (CMRelementToRowIndex(secondMarker1) == CMRelementToRowIndex(secondMarker2)))
  {
    return CMR_ERROR_INPUT;
  }
  if (CMRelementIsColumn(secondMarker1) && CMRelementIsColumn(secondMarker2)
    && (CMRelementToColumnIndex(secondMarker1) == CMRelementToColumnIndex(secondMarker2)))
  {
    return CMR_ERROR_INPUT;
  }

  /* If the first markers are column,row then we exchange both pairs such that the first markers are row,column. */
  if (CMRelementIsColumn(firstMarker1) && CMRelementIsRow(firstMarker2))
  {
    CMR_ELEMENT temp = firstMarker1;
    firstMarker1 = firstMarker2;
    firstMarker2 = temp;
    temp = secondMarker1;
    secondMarker1 = secondMarker2;
    secondMarker2 = temp;
  }

  CMRdbgMsg(0, "Computing 3-sum of %zux%zu and %zux%zu matrix with markers %s,%s for first and %s,%s for second.\n",
    first->numRows, first->numColumns, second->numRows, second->numColumns,
    CMRelementIsRow(firstMarker1) ? "row" : "column", CMRelementIsRow(firstMarker2) ? "row" : "column",
    CMRelementIsRow(secondMarker1) ? "row" : "column", CMRelementIsRow(secondMarker2) ? "row" : "column");

  bool bottomLeft = CMRelementIsRow(firstMarker1) || CMRelementIsRow(firstMarker2);
  bool topRight = CMRelementIsColumn(firstMarker1) || CMRelementIsColumn(firstMarker2);
  CMRdbgMsg(2, "Non-zero bottom-left: %s, non-zero top-right: %s\n", bottomLeft ? "yes" : "no", topRight ? "yes" : "no");

  size_t firstRowMarkers[2] = { SIZE_MAX, SIZE_MAX }; /*< row indices of A where to find a_i (or SIZE_MAX). */
  size_t secondRowMarkers[2] = { SIZE_MAX, SIZE_MAX }; /*< row indices of B where to find b_i (or SIZE_MAX). */
  size_t firstColumnMarkers[2] = { SIZE_MAX, SIZE_MAX }; /*< column indices of A where to find a_i (or SIZE_MAX). */
  size_t secondColumnMarkers[2] = { SIZE_MAX, SIZE_MAX }; /*< column indices of B where to find b_i (or SIZE_MAX). */
  char* markerColumn[2] = { NULL, NULL }; /* Nonzero entries of the column vector among a_i,b_i. */
  size_t markerColumnNumNonzeros[2] = { 0, 0 }; /* Number of nonzeros of markerColumn[i]. */
  size_t markerRowNumNonzeros[2] = { 0, 0 }; /* Number of nonzeros of a_i that is not markerColumn[i]. */

  /* Step 1: fill marker columns. */

  /* Step 1a: rank(bottom-left)-many columns of B are stored in markerColumn. */
  uint8_t bottomLeftRank = (!topRight) ? 2 : (bottomLeft ? 1 : 0);
  for (uint8_t i = 0; i < bottomLeftRank; ++i)
  {
    firstRowMarkers[i] = CMRelementToRowIndex(i == 0 ? firstMarker1 : firstMarker2);
    markerRowNumNonzeros[i] = first->rowSlice[firstRowMarkers[i] + 1] - first->rowSlice[firstRowMarkers[i]];
    CMRdbgMsg(2, "firstRowMarkers[%d] = %zu\n", i, firstRowMarkers[i]);

    secondColumnMarkers[i] = CMRelementToColumnIndex(i == 0 ? secondMarker1 : secondMarker2);
    CMR_CALL( CMRallocStackArray(cmr, &markerColumn[i], second->numColumns) );
    for (size_t row = 0; row < second->numRows; ++row)
    {
      size_t entry;
      CMR_CALL( CMRchrmatFindEntry(second, row, secondColumnMarkers[i], &entry) );
      if (entry == SIZE_MAX)
      {
        markerColumn[i][row] = 0;
        CMRdbgMsg(2, "markerColumn[%d][%zu] is 0.\n", i, row);
      }
      else
      {
        assert(second->entryColumns[entry] == secondColumnMarkers[i]);
        markerColumn[i][row] = second->entryValues[entry];
        CMRdbgMsg(2, "markerColumn[%d][%zu] = %d.\n", i, row, second->entryValues[entry]);
        ++markerColumnNumNonzeros[i];
      }
    }
  }

  /* Step 1b: rank(top-right)-many columns of A are stored in markerColumn. */
  uint8_t topRightRank = 2 - bottomLeftRank;
  for (uint8_t j = 0; j < topRightRank; ++j)
  {
    uint8_t i = bottomLeftRank + j;
    secondRowMarkers[i] = CMRelementToRowIndex(i == 0 ? secondMarker1 : secondMarker2);
    markerRowNumNonzeros[i] = second->rowSlice[secondRowMarkers[i] + 1] - second->rowSlice[secondRowMarkers[i]];

    firstColumnMarkers[i] = CMRelementToColumnIndex(i == 0 ? firstMarker1 : firstMarker2);
    CMR_CALL( CMRallocStackArray(cmr, &markerColumn[i], second->numColumns) );
    for (size_t row = 0; row < first->numRows; ++row)
    {
      size_t entry;
      CMR_CALL( CMRchrmatFindEntry(first, row, firstColumnMarkers[i], &entry) );
      if (entry == SIZE_MAX)
      {
        markerColumn[i][row] = 0;
        CMRdbgMsg(2, "markerColumn[%d][%ld] is 0.\n", i, row);
      }
      else
      {
        assert(first->entryColumns[entry] == firstColumnMarkers[i]);
        markerColumn[i][row] = first->entryValues[entry];
        CMRdbgMsg(2, "markerColumn[%d][%ld] = %d.\n", i, row, first->entryValues[entry]);
        ++markerColumnNumNonzeros[i];
      }
    }
  }

  assert(markerColumn[0]);
  assert(markerColumn[1]);

  /* Step 2: create mapping from old rows/columns to new ones. */
  size_t* firstRowToResultRow = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &firstRowToResultRow, first->numRows) );
  size_t resultRow = 0;
  for (size_t row = 0; row < first->numRows; ++row)
  {
    if ((row == firstRowMarkers[0]) || (row == firstRowMarkers[1]))
      firstRowToResultRow[row] = SIZE_MAX;
    else
      firstRowToResultRow[row] = resultRow++;
  }
  size_t* secondRowToResultRow = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &secondRowToResultRow, second->numRows) );
  for (size_t row = 0; row < second->numRows; ++row)
  {
    if ((row == secondRowMarkers[0]) || (row == secondRowMarkers[1]))
      secondRowToResultRow[row] = SIZE_MAX;
    else
      secondRowToResultRow[row] = resultRow++;
  }
  CMRdbgMsg(2, "resultRow = %zu\n", resultRow);
  assert(resultRow == first->numRows + second->numRows - 2);

  size_t* firstColumnToResultColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &firstColumnToResultColumn, first->numColumns) );
  size_t resultColumn = 0;
  for (size_t column = 0; column < first->numColumns; ++column)
  {
    if ((column == firstColumnMarkers[0]) || (column == firstColumnMarkers[1]))
      firstColumnToResultColumn[column] = SIZE_MAX;
    else
      firstColumnToResultColumn[column] = resultColumn++;
  }
  size_t* secondColumnToResultColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &secondColumnToResultColumn, second->numColumns) );
  for (size_t column = 0; column < second->numColumns; ++column)
  {
    if ((column == secondColumnMarkers[0]) || (column == secondColumnMarkers[1]))
      secondColumnToResultColumn[column] = SIZE_MAX;
    else
      secondColumnToResultColumn[column] = resultColumn++;
  }

  /* Step 3: create the actual matrix. */
  CMRdbgMsg(0, "upper bound on #nonzeros is %zu + %zu + (%zu * %zu - 1) + (%zu * %zu - 1)\n", first->numNonzeros,
    second->numNonzeros, markerRowNumNonzeros[0] - 1, markerColumnNumNonzeros[0] - 1, markerRowNumNonzeros[1] - 1,
    markerColumnNumNonzeros[1] - 1);
  CMR_CALL( CMRchrmatCreate(cmr, presult, first->numRows + second->numRows - 2,
    first->numColumns + second->numColumns - 2,
    first->numNonzeros + second->numNonzeros + (markerRowNumNonzeros[0] - 1) * (markerColumnNumNonzeros[0] - 1)
    + (markerRowNumNonzeros[1] - 1) * (markerColumnNumNonzeros[1] - 1) - 1) );
  CMR_CHRMAT* result = *presult;
  resultRow = 0;
  size_t resultNonzero = 0;
  CMR_ERROR error = CMR_OKAY;

  /* Step 3a: top part. */
  for (size_t row = 0; row < first->numRows; ++row)
  {
    if ((row == firstRowMarkers[0]) || (row == firstRowMarkers[1]))
      continue;

    result->rowSlice[resultRow] = resultNonzero;
    CMRdbgMsg(2, "Row %zu starts at nonzero #%zu and belongs to top part.\n", resultRow, resultNonzero);

    /* First matrix is top-left. */
    size_t begin = first->rowSlice[row];
    size_t beyond = first->rowSlice[row + 1];
    for (size_t e = begin; e < beyond; ++e)
    {
      size_t column = first->entryColumns[e];
      if (column == firstColumnMarkers[0] || column == firstColumnMarkers[1])
        continue;
      result->entryValues[resultNonzero] = first->entryValues[e];
      result->entryColumns[resultNonzero] = firstColumnToResultColumn[column];
      CMRdbgMsg(4, "Adding nonzero #%zu at %zu,%zu with value %d from top-left.\n", resultNonzero, resultRow,
        result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
      ++resultNonzero;
    }

    if (topRight && bottomLeft)
    {
      /* Compute values of rank-1 top-right matrix. */
      if (markerColumn[1][row] != 0)
      {
        size_t begin = second->rowSlice[secondRowMarkers[1]];
        size_t beyond = second->rowSlice[secondRowMarkers[1] + 1];
        for (size_t e = begin; e < beyond; ++e)
        {
          if (second->entryColumns[e] == secondColumnMarkers[0])
            continue;

          /* Take care of modulo operations. */
          int resultValue = (int)markerColumn[1][row] * (int)second->entryValues[e];
          assert(resultValue != 0);
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }
          assert(resultValue != 0);

          if (resultValue > INT8_MAX || resultValue < INT8_MIN)
          {
            error = CMR_ERROR_OVERFLOW;
            goto cleanup;
          }

          result->entryValues[resultNonzero] = resultValue;
          result->entryColumns[resultNonzero] = secondColumnToResultColumn[second->entryColumns[e]];
          CMRdbgMsg(4, "Adding nonzero #%zu at %zu,%zu with value %d from product in top-right.\n", resultNonzero, resultRow,
            result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
          ++resultNonzero;
        }
      }
    }
    else if (topRight && !bottomLeft)
    {
      /* Compute values of rank-2 top-right matrix. */
      if (markerColumn[0][row] != 0 || markerColumn[1][row] != 0)
      {
        size_t e1 = second->rowSlice[secondRowMarkers[0]];
        size_t beyond1 = second->rowSlice[secondRowMarkers[0] + 1];
        size_t column1 = (e1 < beyond1) ? second->entryColumns[e1] : SIZE_MAX;
        size_t e2 = second->rowSlice[secondRowMarkers[1]];
        size_t beyond2 = second->rowSlice[secondRowMarkers[1] + 1];
        size_t column2 = (e2 < beyond2) ? second->entryColumns[e2] : SIZE_MAX;
        while (column1 != SIZE_MAX || column2 != SIZE_MAX)
        {
          int resultValue;
          size_t nextColumn;
          if (column1 < column2)
          {
            resultValue = (int)markerColumn[0][row] * (int)second->entryValues[e1];
            nextColumn = column1;
            ++e1;
            column1 = (e1 < beyond1) ? second->entryColumns[e1] : SIZE_MAX;
          }
          else if (column2 < column1)
          {
            resultValue = (int)markerColumn[1][row] * (int)second->entryValues[e2];
            nextColumn = column2;
            ++e2;
            column2 = (e2 < beyond2) ? second->entryColumns[e2] : SIZE_MAX;
          }
          else
          {
            resultValue = (int)markerColumn[0][row] * (int)second->entryValues[e1]
              + (int)markerColumn[1][row] * (int)second->entryValues[e2];
            nextColumn = column1;
            ++e1;
            column1 = (e1 < beyond1) ? second->entryColumns[e1] : SIZE_MAX;
            ++e2;
            column2 = (e2 < beyond2) ? second->entryColumns[e2] : SIZE_MAX;
          }

          /* Take care of modulo operations. */
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }

          if (resultValue > INT8_MAX || resultValue < INT8_MIN)
          {
            error = CMR_ERROR_OVERFLOW;
            goto cleanup;
          }

          if (resultValue != 0)
          {
            result->entryValues[resultNonzero] = resultValue;
            result->entryColumns[resultNonzero] = secondColumnToResultColumn[nextColumn];
            CMRdbgMsg(4, "Adding nonzero #%zu at %zu,%zu with value %d.\n", resultNonzero, resultRow,
              result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
            ++resultNonzero;
          }
          else
          {
            CMRdbgMsg(4, "Computed a zero at %zu,%zu in the rank-2 block.\n", resultRow, resultColumn);
          }
        }
      }
    }

    ++resultRow;
  }

  /* Step 3b: bottom part. */
  for (size_t row = 0; row < second->numRows; ++row)
  {
    if ((row == secondRowMarkers[0]) || (row == secondRowMarkers[1]))
      continue;

    result->rowSlice[resultRow] = resultNonzero;
    CMRdbgMsg(2, "Row %zu starts at nonzero #%zu.\n", resultRow, resultNonzero);

    if (bottomLeft && topRight)
    {
      /* Compute values of rank-1 bottom-left matrix. */
      if (markerColumn[0][row] != 0)
      {
        size_t begin = first->rowSlice[firstRowMarkers[0]];
        size_t beyond = first->rowSlice[firstRowMarkers[0] + 1];
        for (size_t e = begin; e < beyond; ++e)
        {
          if (first->entryColumns[e] == firstColumnMarkers[1])
            continue;

          int resultValue = (int)markerColumn[0][row] * (int)first->entryValues[e];
          assert(resultValue != 0);

          /* Take care of modulo operations. */
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }
          assert(resultValue != 0);

          if (resultValue > INT8_MAX || resultValue < INT8_MIN)
          {
            error = CMR_ERROR_OVERFLOW;
            goto cleanup;
          }

          result->entryValues[resultNonzero] = resultValue;
          result->entryColumns[resultNonzero] = firstColumnToResultColumn[first->entryColumns[e]];
          CMRdbgMsg(4, "Adding nonzero #%zu at %zu,%zu with value %d from bottom-left.\n", resultNonzero, resultRow,
            result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
          ++resultNonzero;
        }
      }
    }
    else if (!topRight)
    {
      /* Compute values of rank-2 bottom-left matrix. */
      if (markerColumn[0][row] != 0 || markerColumn[1][row] != 0)
      {
        size_t e1 = first->rowSlice[firstRowMarkers[0]];
        size_t beyond1 = first->rowSlice[firstRowMarkers[0] + 1];
        size_t column1 = (e1 < beyond1) ? first->entryColumns[e1] : SIZE_MAX;
        size_t e2 = first->rowSlice[firstRowMarkers[1]];
        size_t beyond2 = first->rowSlice[firstRowMarkers[1] + 1];
        size_t column2 = (e2 < beyond2) ? first->entryColumns[e2] : SIZE_MAX;
        while (column1 != SIZE_MAX || column2 != SIZE_MAX)
        {
          int resultValue;
          size_t nextColumn;
          if (column1 < column2)
          {
            resultValue = (int)markerColumn[0][row] * (int)first->entryValues[e1];
            CMRdbgMsg(4, "Calculating entry %zu,%zu (was %zu,%zu) as product of %d and %d = %d.\n", resultRow, row,
              firstColumnToResultColumn[column1], column1, markerColumn[0][row], first->entryValues[e1], resultValue);
            nextColumn = column1;
            ++e1;
            column1 = (e1 < beyond1) ? first->entryColumns[e1] : SIZE_MAX;
          }
          else if (column2 < column1)
          {
            resultValue = (int)markerColumn[1][row] * (int)first->entryValues[e2];
            nextColumn = column2;
            ++e2;
            column2 = (e2 < beyond2) ? first->entryColumns[e2] : SIZE_MAX;
          }
          else
          {
            resultValue = (int)markerColumn[0][row] * (int)first->entryValues[e1]
              + (int)markerColumn[1][row] * (int)first->entryValues[e2];
            nextColumn = column1;
            ++e1;
            column1 = (e1 < beyond1) ? first->entryColumns[e1] : SIZE_MAX;
            ++e2;
            column2 = (e2 < beyond2) ? first->entryColumns[e2] : SIZE_MAX;
          }

          /* Take care of modulo operations. */
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }

          if (resultValue > INT8_MAX || resultValue < INT8_MIN)
          {
            error = CMR_ERROR_OVERFLOW;
            goto cleanup;
          }

          if (resultValue != 0)
          {
            result->entryValues[resultNonzero] = resultValue;
            result->entryColumns[resultNonzero] = firstColumnToResultColumn[nextColumn];
            CMRdbgMsg(4, "Adding nonzero #%zu at %ld,%ld with value %d from bottom-left.\n", resultNonzero, resultRow,
              result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
            ++resultNonzero;
          }
          else
          {
            CMRdbgMsg(4, "Computed a zero at %ld,%ld in the rank-2 block.\n", resultRow, resultColumn);
          }
        }
      }
    }

    /* Second matrix is bottom-right. */
    for (size_t e = second->rowSlice[row]; e < second->rowSlice[row + 1]; ++e)
    {
      size_t column = second->entryColumns[e];
      if ((column == secondColumnMarkers[0]) || (column == secondColumnMarkers[1]))
        continue;
      result->entryValues[resultNonzero] = second->entryValues[e];
      result->entryColumns[resultNonzero] = secondColumnToResultColumn[column];
      CMRdbgMsg(4, "Adding nonzero #%zu at %ld,%ld with value %d from bottom-right.\n", resultNonzero, resultRow,
        result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
      ++resultNonzero;
    }

    ++resultRow;
  }

  CMRdbgMsg(2, "Allocated memory for %zu nonzeros of which %zu are used.\n", result->numNonzeros, resultNonzero);

  result->rowSlice[result->numRows] = resultNonzero;
  assert(result->numNonzeros >= resultNonzero);
  result->numNonzeros = resultNonzero;

  CMRconsistencyAssert( CMRchrmatConsistency(result) );

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &secondColumnToResultColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &firstColumnToResultColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &secondRowToResultRow) );
  CMR_CALL( CMRfreeStackArray(cmr, &firstRowToResultRow) );
  CMR_CALL( CMRfreeStackArray(cmr, &markerColumn[1]) );
  CMR_CALL( CMRfreeStackArray(cmr, &markerColumn[0]) );

  return error;
}

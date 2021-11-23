// #define CMR_DEBUG

#include <cmr/ctu.h>

#include <cmr/tu.h>

#include "env_internal.h"

#include <assert.h>
#include <stdint.h>

CMR_ERROR CMRcomplementRowColumn(CMR* cmr, CMR_CHRMAT* matrix, size_t complementRow, size_t complementColumn,
  CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(complementRow < matrix->numRows || complementRow == SIZE_MAX);
  assert(complementColumn < matrix->numColumns || complementColumn == SIZE_MAX);
  assert(presult);

  /* We copy directly if no complementing is requested. */

  if (complementRow == SIZE_MAX && complementColumn == SIZE_MAX)
  {
    CMR_CALL( CMRchrmatCopy(cmr, matrix, presult) );

    return CMR_OKAY;
  }

  /* Create a dense copy of the column to be complemented. */
  char* complementColumnEntries = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &complementColumnEntries, matrix->numRows) );
  if (complementColumn < SIZE_MAX)
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
    {
      size_t e;
      CMR_CALL( CMRchrmatFindEntry(matrix, row, complementColumn, &e) );
      complementColumnEntries[row] = (e == SIZE_MAX) ? 0 : 1;
    }
  }
  else
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
      complementColumnEntries[row] = 0;
  }

  /* Create a dense copy of the row to be complemented. */
  char* complementRowEntries = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &complementRowEntries, matrix->numColumns) );
  if (complementRow < SIZE_MAX)
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
      complementRowEntries[column] = 0;
    size_t first = matrix->rowSlice[complementRow];
    size_t beyond = matrix->rowSlice[complementRow + 1];
    for (size_t entry = first; entry < beyond; ++entry)
      complementRowEntries[matrix->entryColumns[entry]] = 1;
  }
  else
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
      complementRowEntries[column] = 0;
  }

  /* Indicator for entry in complementRow, complementColumn. */
  char complementRowColumn1 = complementColumn < SIZE_MAX ? complementRowEntries[complementColumn] : 0;

  /* Swipe over the matrix to create the complemented one. */
  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, 0) );
  CMR_CHRMAT* result = *presult;
  size_t memNonzeros = matrix->numNonzeros + 256;
  CMR_CALL( CMRallocBlockArray(cmr, &result->entryColumns, memNonzeros) );
  size_t resultEntry = 0;
  result->rowSlice[0] = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t matrixEntry = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    size_t matrixColumn = (matrixEntry < beyond) ? matrix->entryColumns[matrixEntry] : SIZE_MAX;
    for (size_t column = 0; column < matrix->numColumns; ++column)
    {
      bool isNonzero = (column == matrixColumn);
      if (isNonzero)
      {
        ++matrixEntry;
        matrixColumn = (matrixEntry < beyond) ? matrix->entryColumns[matrixEntry] : SIZE_MAX;
      }

      if (row == complementRow)
      {
        if (column != complementColumn && complementRowColumn1)
          isNonzero = !isNonzero;
      }
      else
      {
        if (column == complementColumn)
        {
          if (complementRowColumn1)
            isNonzero = !isNonzero;
        }
        else
        {
          if (complementRowColumn1 + complementColumnEntries[row] + complementRowEntries[column] % 2 == 1)
            isNonzero = !isNonzero;
        }
      }

      if (isNonzero)
      {
        if (resultEntry == memNonzeros)
        {
          memNonzeros *= 2;
          CMR_CALL( CMRreallocBlockArray(cmr, &result->entryColumns, memNonzeros) );
        }
        result->entryColumns[resultEntry++] = column;
      }
    }
    result->rowSlice[row + 1] = resultEntry;
  }

  /* Finalize the size of entryColumns, allocate entryValues for the first time and fill the latter with 1s. */

  if (resultEntry < memNonzeros)
    CMR_CALL( CMRchrmatChangeNumNonzeros(cmr, result, resultEntry) );
  for (size_t e = 0; e < resultEntry; ++e)
    result->entryValues[e] = 1;

  CMR_CALL( CMRfreeStackArray(cmr, &complementRowEntries) );
  CMR_CALL( CMRfreeStackArray(cmr, &complementColumnEntries) );

  return CMR_OKAY;
}

CMR_ERROR CMRtestComplementTotalUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisComplementTotallyUnimodular,
  size_t* pcomplementRow, size_t* pcomplementColumn)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(pisComplementTotallyUnimodular);

  /* Create a dense version on the stack. */

  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;
  char* dense = NULL;
  size_t sizeDense = numRows * numColumns;
  CMR_CALL( CMRallocStackArray(cmr, &dense, sizeDense) );
  for (size_t i = 0; i < sizeDense; ++i)
    dense[i] = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
      dense[numColumns * row + matrix->entryColumns[entry]] = matrix->entryValues[entry];
  }

  /* Create complemented matrices and call test for total unimodularity. */
  CMR_CHRMAT* complementedMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &complementedMatrix, numRows, numColumns, sizeDense) );
  for (size_t i = 0; i < sizeDense; ++i)
    complementedMatrix->entryValues[i] = 1;
  *pisComplementTotallyUnimodular = true;
  for (size_t complementRow = 0; complementRow <= numRows && *pisComplementTotallyUnimodular; ++complementRow)
  {
    bool hasComplementRow = complementRow < numRows;
    for (size_t complementColumn = 0; complementColumn <= numColumns; ++complementColumn)
    {
      bool hasComplementColumn = complementColumn < numColumns;

      /* Indicator for entry in complementRow, complementColumn. */
      char complementRowColumn1 = hasComplementRow && hasComplementColumn ?
        dense[numColumns * complementRow + complementColumn] : 0;

      /* Fill complemented matrix as a sparse char matrix. */
      complementedMatrix->numNonzeros = 0;
      complementedMatrix->rowSlice[0] = 0;
      for (size_t row = 0; row < numRows; ++row)
      {
        for (size_t column = 0; column < numColumns; ++column)
        {
          bool isNonzero = dense[numColumns * row + column];
          if (row == complementRow)
          {
            if (column != complementColumn && complementRowColumn1)
              isNonzero = !isNonzero;
          }
          else
          {
            if (column == complementColumn)
            {
              if (complementRowColumn1)
                isNonzero = !isNonzero;
            }
            else
            {
              if ((complementRowColumn1 + (hasComplementColumn ? dense[numColumns * row + complementColumn] : 0)
                + (hasComplementRow ? dense[numColumns * complementRow + column] : 0)) % 2 == 1)
              {
                isNonzero = !isNonzero;
              }
            }
          }

          if (isNonzero)
          {
            complementedMatrix->entryColumns[complementedMatrix->numNonzeros] = column;
            complementedMatrix->numNonzeros++;
          }
        }
        complementedMatrix->rowSlice[row + 1] = complementedMatrix->numNonzeros;
      }

      bool isTU = false;
      CMR_CALL( CMRtestTotalUnimodularity(cmr, complementedMatrix, &isTU, NULL, NULL, NULL) );
      if (!isTU)
      {
        if (pcomplementRow)
          *pcomplementRow = complementRow;
        if (pcomplementColumn)
          *pcomplementColumn = complementColumn;

        *pisComplementTotallyUnimodular = false;
        break;
      }
    }
  }
  CMR_CALL( CMRchrmatFree(cmr, &complementedMatrix) );

  CMR_CALL( CMRfreeStackArray(cmr, &dense) );

  return CMR_OKAY;
}

// #define CMR_DEBUG /* Uncomment to debug. */

#include "hereditary_property.h"
#include <stdint.h>
#include <time.h>

CMR_ERROR CMRtestHereditaryPropertySimple(CMR* cmr, CMR_CHRMAT* matrix, HereditaryPropertyTest testFunction,
  void* testData, CMR_SUBMAT** psubmatrix, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(testFunction);
  assert(psubmatrix);

  clock_t time = clock();
  size_t* essentialRows = NULL;
  size_t numEssentialRows = 0;
  CMR_CALL( CMRallocStackArray(cmr, &essentialRows, matrix->numRows) );
  size_t* essentialColumns = NULL;
  size_t numEssentialColumns = 0;
  CMR_CALL( CMRallocStackArray(cmr, &essentialColumns, matrix->numColumns) );
  
  CMR_ELEMENT* candidates = NULL;
  size_t numCandidates = matrix->numRows + matrix->numColumns;
  CMR_CALL( CMRallocStackArray(cmr, &candidates, numCandidates) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    candidates[row] = CMRrowToElement(row);
  for (size_t column = 0; column < matrix->numColumns; ++column)
    candidates[matrix->numRows + column] = CMRcolumnToElement(column);

  CMR_CHRMAT* current = NULL;
  CMR_CALL( CMRchrmatCopy(cmr, matrix, &current) );

  CMR_CHRMAT* candidateMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &candidateMatrix, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  candidateMatrix->numNonzeros = 0;

  while (numCandidates > 0)
  {
    CMR_ELEMENT candidateElement = candidates[--numCandidates];
    size_t removedRow = SIZE_MAX;
    size_t removedColumn = SIZE_MAX;
    if (CMRelementIsRow(candidateElement))
      removedRow = CMRelementToRowIndex(candidateElement);
    else
      removedColumn = CMRelementToColumnIndex(candidateElement);

    /* Copy current matrix to candidate matrix, except for removed row/column. */
    candidateMatrix->numNonzeros = 0;
    for (size_t row = 0; row < matrix->numRows; ++row)
    {
      candidateMatrix->rowSlice[row] = candidateMatrix->numNonzeros;
      if (row == removedRow)
        continue;
      size_t first = current->rowSlice[row];
      size_t beyond = current->rowSlice[row + 1];
      for (size_t e = first; e < beyond; ++e)
      {
        size_t column = current->entryColumns[e];
        if (column == removedColumn)
          continue;
        candidateMatrix->entryColumns[candidateMatrix->numNonzeros] = column;
        candidateMatrix->entryValues[candidateMatrix->numNonzeros] = current->entryValues[e];
        candidateMatrix->numNonzeros++;
      }
    }
    candidateMatrix->rowSlice[matrix->numRows] = candidateMatrix->numNonzeros;

    /* Invoke test. */
    bool hasProperty;
    CMR_SUBMAT* submatrix = NULL;
    double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    if (remainingTime < 0)
    {
      CMR_CALL( CMRchrmatFree(cmr, &candidateMatrix) );
      CMR_CALL( CMRchrmatFree(cmr, &current) );
      CMR_CALL( CMRfreeStackArray(cmr, &candidates) );
      CMR_CALL( CMRfreeStackArray(cmr, &essentialColumns) );
      CMR_CALL( CMRfreeStackArray(cmr, &essentialRows) );
      return CMR_ERROR_TIMEOUT;
    }

    CMRdbgMsg(2, "\n!!! Hereditary property test queries the test oracle!!!\n\n");
    CMR_ERROR error = testFunction(cmr, candidateMatrix, testData, &hasProperty, &submatrix, remainingTime);
    if (error == CMR_ERROR_TIMEOUT)
    {
      CMR_CALL( CMRchrmatFree(cmr, &candidateMatrix) );
      CMR_CALL( CMRchrmatFree(cmr, &current) );
      CMR_CALL( CMRfreeStackArray(cmr, &candidates) );
      CMR_CALL( CMRfreeStackArray(cmr, &essentialColumns) );
      CMR_CALL( CMRfreeStackArray(cmr, &essentialRows) );
      return CMR_ERROR_TIMEOUT;
    }
    CMR_CALL(error);

    CMRdbgMsg(2, "\n!!! Property %s present.\n\n", hasProperty ? "IS" : "is NOT");

    assert(!submatrix); // TODO: we cannot deal with this, yet.

    if (hasProperty)
    {
      if (removedRow < SIZE_MAX)
        essentialRows[numEssentialRows++] = removedRow;
      else
        essentialColumns[numEssentialColumns++] = removedColumn;
    }
    else
    {
      /* Swap candidateMatrix and current. */
      CMR_CHRMAT* temp = candidateMatrix;
      candidateMatrix = current;
      current = temp;
    }
  }

  CMR_CALL( CMRchrmatFree(cmr, &candidateMatrix) );
  CMR_CALL( CMRchrmatFree(cmr, &current) );

  /* Extract the submatrix. */
  CMR_CALL( CMRsubmatCreate(cmr, numEssentialRows, numEssentialColumns, psubmatrix) );
  CMR_SUBMAT* submatrix = *psubmatrix;
  for (size_t row = 0; row < submatrix->numRows; ++row)
    submatrix->rows[row] = essentialRows[row];
  for (size_t column = 0; column < submatrix->numColumns; ++column)
    submatrix->columns[column] = essentialColumns[column];

  CMR_CALL( CMRfreeStackArray(cmr, &candidates) );
  CMR_CALL( CMRfreeStackArray(cmr, &essentialColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &essentialRows) );

  return CMR_OKAY;
}

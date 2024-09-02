// #define CMR_DEBUG /* Uncomment to debug. */
// #define CMR_DEBUG_MATRICES /* Uncomment to debug the actual matrices. */

#include "hereditary_property.h"
#include <stdint.h>
#include <time.h>

/**
 * \brief Tests a submatrix of the given \p current matrix for the hereditary property defined by a given
 *        \p testFunction.
 *
 * The algorithm finds the submatrix by successively single zeroing out rows or columns. The submatrix is specified by
 * \p remainingRows and \p remainingColumns. The \p current matrix is modified and later deleted.
 */

static
CMR_ERROR testHereditaryPropertyNaive(
  CMR* cmr,                             /**< \ref CMR environment. */
  CMR_CHRMAT* current,                  /**< Some matrix not having the hereditary property. */
  HereditaryPropertyTest testFunction,  /**< Test function. */
  void* testData,                       /**< Data to be forwarded to the test function. */
  CMR_SUBMAT** psubmatrix,              /**< Pointer for storing a minimal submatrix not having the property. */
  size_t numRemainingRows,              /**< Number of rows of submatrix to search in. */
  size_t* remainingRows,                /**< Rows of submatrix to search in. */
  size_t numRemainingColumns,           /**< Number of columns of submatrix to search in. */
  size_t* remainingColumns,             /**< Columns of submatrix to search in. */
  double timeLimit                      /**< Time limit to impose. */
  )
{
  assert(cmr);
  assert(current);
  assert(testFunction);
  assert(numRemainingRows > 0);
  assert(remainingRows || numRemainingRows == current->numRows);
  assert(numRemainingColumns > 0);
  assert(remainingColumns || numRemainingColumns == current->numColumns);

  CMR_ERROR error = CMR_OKAY;
  clock_t startTime = clock();
  double remainingTime;

  CMR_CHRMAT* candidateMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &candidateMatrix, current->numRows, current->numColumns, current->numNonzeros) );
  candidateMatrix->numNonzeros = 0;

  size_t* essentialRows = NULL;
  size_t numEssentialRows = 0;
  CMR_CALL( CMRallocStackArray(cmr, &essentialRows, numRemainingRows) );
  size_t* essentialColumns = NULL;
  size_t numEssentialColumns = 0;
  CMR_CALL( CMRallocStackArray(cmr, &essentialColumns, numRemainingColumns) );

  CMRdbgMsg(0, "Simply searching minimal violator of %zux%zu matrix.\n", numRemainingRows, numRemainingColumns);

  size_t numRemainingElements = numRemainingRows + numRemainingColumns;
  for (size_t i = numRemainingElements; i > 0; --i)
  {
    size_t removedElement = i - 1;
    size_t removedRow = SIZE_MAX;
    size_t removedColumn = SIZE_MAX;
    if (removedElement < numRemainingColumns)
      removedColumn = remainingColumns ? remainingColumns[removedElement] : removedElement;
    else
    {
      removedRow = remainingRows ? remainingRows[removedElement - numRemainingColumns]
        : (removedElement - numRemainingColumns);
    }

    if (removedElement < numRemainingColumns)
      CMRdbgMsg(2, "Removing c%zu.\n", removedColumn + 1);
    else
      CMRdbgMsg(2, "Removing r%zu.\n", removedRow + 1);

    /* Copy current matrix to candidate matrix, except for removed row/column. */
    candidateMatrix->numNonzeros = 0;
    for (size_t row = 0; row < current->numRows; ++row)
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
    candidateMatrix->rowSlice[current->numRows] = candidateMatrix->numNonzeros;

    CMRdbgMsg(4, "Submatrix has %zu nonzeros.\n", candidateMatrix->numNonzeros);

    /* Invoke test. */
    bool hasProperty;
    CMR_SUBMAT* submatrix = NULL;
    remainingTime = timeLimit - (clock() - startTime) * 1.0 / CLOCKS_PER_SEC;
    if (remainingTime < 0)
    {
      error = CMR_ERROR_TIMEOUT;
      goto cleanup;
    }

    CMRdbgMsg(4, "\n!!!     Hereditary property test queries the test oracle!!!\n\n");
    error = testFunction(cmr, candidateMatrix, testData, &hasProperty, &submatrix, remainingTime);
    if (error != CMR_OKAY)
      goto cleanup;

    CMRdbgMsg(4, "\n!!! Property %s present.\n\n", hasProperty ? "IS" : "is NOT");

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

#ifdef CMR_DEBUG_MATRICES

  CMRdbgMsg(2, "Final zeroed matrix:\n");
  CMRchrmatPrintDense(cmr, current, stdout, '0', true);

#endif /* CMR_DEBUG_MATRICES */

  CMRdbgMsg(2, "Extracting a %zux%zu submatrix.\n", numEssentialRows, numEssentialColumns);

  /* Extract the submatrix. */
  CMR_CALL( CMRsubmatCreate(cmr, numEssentialRows, numEssentialColumns, psubmatrix) );
  CMR_SUBMAT* submatrix = *psubmatrix;
  for (size_t row = 0; row < submatrix->numRows; ++row)
    submatrix->rows[row] = essentialRows[row];
  for (size_t column = 0; column < submatrix->numColumns; ++column)
    submatrix->columns[column] = essentialColumns[column];

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &essentialColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &essentialRows) );
  CMR_CALL( CMRchrmatFree(cmr, &candidateMatrix) );
  CMR_CALL( CMRchrmatFree(cmr, &current) );

  return error;
}

CMR_ERROR CMRtestHereditaryPropertyNaive(CMR* cmr, CMR_CHRMAT* matrix, HereditaryPropertyTest testFunction,
  void* testData, CMR_SUBMAT** psubmatrix, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(testFunction);
  assert(psubmatrix);

  CMR_CHRMAT* current = NULL;
  CMR_CALL( CMRchrmatCopy(cmr, matrix, &current) );

  CMR_ERROR error = testHereditaryPropertyNaive(cmr, current, testFunction, testData, psubmatrix, current->numRows,
    NULL, current->numColumns, NULL, timeLimit);

  return error;
}

CMR_ERROR CMRtestHereditaryPropertyGreedy(CMR* cmr, CMR_CHRMAT* matrix, HereditaryPropertyTest testFunction,
  void* testData, CMR_SUBMAT** psubmatrix, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(testFunction);
  assert(psubmatrix);

  clock_t startTime = clock();
  double remainingTime;

  /* Remaining rows/column are those that we have restricted ourself to already. */
  size_t* remainingRows = NULL;
  size_t numRemainingRows = matrix->numRows;
  CMR_CALL( CMRallocStackArray(cmr, &remainingRows, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    remainingRows[row] = row;
  size_t* remainingColumns = NULL;
  size_t numRemainingColumns = matrix->numColumns;
  CMR_CALL( CMRallocStackArray(cmr, &remainingColumns, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    remainingColumns[column] = column;
  CMR_CHRMAT* remainingMatrix = NULL;
  CMR_CALL( CMRchrmatCopy(cmr, matrix, &remainingMatrix) );

  /* Candidate rows/columns are a subset of the remaining ones that we inspect. */
  size_t* candidateRows = NULL;
  size_t numCandidateRows;
  CMR_CALL( CMRallocStackArray(cmr, &candidateRows, matrix->numRows) );
  size_t numCandidateColumns;
  bool* columnsCandidate = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnsCandidate, matrix->numColumns) );
  CMR_CHRMAT* candidateMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &candidateMatrix, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  candidateMatrix->numNonzeros = 0;

  CMRdbgMsg(0, "Greedily searching minimal violator of %zux%zu matrix.\n", matrix->numRows, matrix->numColumns);
#ifdef CMR_DEBUG_MATRICES
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
#endif /* CMR_DEBUG_MATRICES */

  CMR_ERROR error = CMR_OKAY;
  for (size_t numChunks = 16; numChunks >= 4; numChunks /= 2)
  {
    CMRdbgMsg(2, "Remaining matrix is %zux%zu.\n", numRemainingRows, numRemainingColumns);

    double fraction = 1.0 / numChunks;
    for (size_t elementChunk = 0; elementChunk < 2 * numChunks; ++elementChunk)
    {
      size_t firstRemainingRow = 0;
      size_t beyondRemainingRow = numRemainingRows;
      size_t firstRemainingColumn = 0;
      size_t beyondRemainingColumn = numRemainingColumns;
      if (elementChunk < numChunks)
      {
        firstRemainingRow = (size_t)(fraction * elementChunk * numRemainingRows);
        if (elementChunk + 1 < numChunks)
          beyondRemainingRow = (size_t)(fraction * (elementChunk + 1) * numRemainingRows);
      }
      else
      {
        firstRemainingColumn = (size_t)(fraction * (elementChunk - numChunks) * numRemainingColumns );
        if (elementChunk + 1 < 2 * numChunks)
          beyondRemainingColumn = (size_t)(fraction * (elementChunk - numChunks + 1) * numRemainingColumns);
      }

      numCandidateRows = beyondRemainingRow - firstRemainingRow;
      numCandidateColumns = beyondRemainingColumn - firstRemainingColumn;

      if (numCandidateRows < 2 || numCandidateColumns < 2)
        continue;

      for (size_t r = 0; r < numCandidateRows; ++r)
        candidateRows[r] = remainingRows[firstRemainingRow + r];
      for (size_t column = 0; column < matrix->numColumns; ++column)
        columnsCandidate[column] = false;
      for (size_t c = firstRemainingColumn; c < beyondRemainingColumn; ++c)
        columnsCandidate[remainingColumns[c]] = true;

#ifdef CMR_DEBUG
      CMRdbgMsg(4, "Chunk #%zu of 2*%zu restricts to %zux%zu matrix.\n", elementChunk, numChunks,
        numCandidateRows, numCandidateColumns);
      CMRdbgMsg(6, "Rows:   ");
      for (size_t r = 0; r < numCandidateRows; ++r)
        printf(" r%zu", candidateRows[r] + 1);
      CMRdbgMsg(0, "\n");
      CMRdbgMsg(6, "Columns:");
      for (size_t column = 0; column < matrix->numColumns; ++column)
      {
        if (columnsCandidate[column])
          printf(" c%zu", column + 1);
      }
      CMRdbgMsg(0, "\n");
#endif /* CMR_DEBUG */

      /* Extract the candidate matrix. */
      candidateMatrix->numNonzeros = 0;
      size_t nextRow = 0;
      for (size_t r = 0; r < numCandidateRows; ++r)
      {
        size_t candidateRow = candidateRows[r];
        for (; nextRow <= candidateRow; ++nextRow)
          candidateMatrix->rowSlice[nextRow] = candidateMatrix->numNonzeros;

        size_t first = remainingMatrix->rowSlice[candidateRow];
        size_t beyond = remainingMatrix->rowSlice[candidateRow + 1];
        for (size_t e = first; e < beyond; ++e)
        {
          size_t column = remainingMatrix->entryColumns[e];
          if (columnsCandidate[column])
          {
            candidateMatrix->entryColumns[candidateMatrix->numNonzeros] = column;
            candidateMatrix->entryValues[candidateMatrix->numNonzeros] = remainingMatrix->entryValues[e];
            candidateMatrix->numNonzeros++;
          }
        }
      }
      for (; nextRow <= matrix->numRows; ++nextRow)
        candidateMatrix->rowSlice[nextRow] = candidateMatrix->numNonzeros;
      CMRdbgMsg(4, "Candidate (sub)matrix has %zu nonzeros.\n", candidateMatrix->numNonzeros);
#ifdef CMR_DEBUG_MATRICES
      CMR_CALL( CMRchrmatPrintDense(cmr, candidateMatrix, stdout, '0', true) );
#endif /* CMR_DEBUG_MATRICES */


      /* Invoke test. */
      bool hasProperty;
      CMR_SUBMAT* submatrix = NULL;
      remainingTime = timeLimit - (clock() - startTime) * 1.0 / CLOCKS_PER_SEC;
      if (remainingTime < 0)
      {
        error = CMR_ERROR_TIMEOUT;
        goto cleanup;
      }

      CMRdbgMsg(4, "\n!!!     Hereditary property test queries the test oracle!!!\n\n");
      error = testFunction(cmr, candidateMatrix, testData, &hasProperty, &submatrix, remainingTime);
      if (error != CMR_OKAY)
        goto cleanup;

      CMRdbgMsg(4, "\n!!! Property %s present.\n\n", hasProperty ? "IS" : "is NOT");

      assert(!submatrix); // TODO: we cannot deal with this, yet.

      if (!hasProperty)
      {
        /* Update remaining rows/columns. */
        numRemainingRows = numCandidateRows;
        for (size_t r = 0; r < numRemainingRows; ++r)
          remainingRows[r] = candidateRows[r];
        numRemainingColumns = numCandidateColumns;
        size_t c = 0;
        for (size_t column = 0; column < matrix->numColumns; ++column)
        {
          if (columnsCandidate[column])
            remainingColumns[c++] = column;
        }

        /* Swap candidateMatrix and remainingMatrix. */
        CMR_CHRMAT* temp = candidateMatrix;
        candidateMatrix = remainingMatrix;
        remainingMatrix = temp;

        /* Ensure that #chunks does not change. */
        numChunks *= 2;
        break;
      }
    }
  }

  remainingTime = timeLimit - (clock() - startTime) * 1.0 / CLOCKS_PER_SEC;
  if (remainingTime < 0)
    error = CMR_ERROR_TIMEOUT;
  else
  {
    error = testHereditaryPropertyNaive(cmr, remainingMatrix, testFunction, testData, psubmatrix, numRemainingRows,
      remainingRows, numRemainingColumns, remainingColumns, remainingTime);
    remainingMatrix = NULL;
  }

cleanup:
  CMR_CALL( CMRchrmatFree(cmr, &candidateMatrix) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnsCandidate) );
  CMR_CALL( CMRfreeStackArray(cmr, &candidateRows) );
  CMR_CALL( CMRchrmatFree(cmr, &remainingMatrix) );
  CMR_CALL( CMRfreeStackArray(cmr, &remainingColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &remainingRows) );

  return error;
}

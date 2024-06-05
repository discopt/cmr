// #define CMR_DEBUG /* Uncomment to debug this file. */

#include "seymour_internal.h"
#include "env_internal.h"

#include "densematrix.h"
#include "hashtable.h"

#include <time.h>

typedef struct
{
  long long hashValue;              /**< \brief Hash value of this element. */
  size_t hashEntry;                 /**< \brief Entry in hashtable. */
  size_t numNonzeros;               /**< \brief Number of nonzeros in (part parallel to) processed submatrix. */
  CMR_ELEMENT representative;       /**< \brief Parallel element. */
  CMR_ELEMENT predecessor;          /**< \brief Predecessor row/column in BFS. */
  bool isProcessed : 1;             /**< \brief Whether this row/column belongs to processed submatrix. */
  bool isSource : 1;                /**< \brief Whether this row/column is a source node in the BFS. */
  bool isTarget : 1;                /**< \brief Whether this row/column is a target node in the BFS. */
  bool isFlipped : 1;               /**< \brief Whether this row/column is a flip node. Entries in a flip row and a
                                                flip column change roles with respect to being edges. */
  bool inQueue : 1;                 /**< \brief Whether this row/column is in the BFS queue. */
} ElementData;

static
CMR_ERROR dbgPrintDenseSequence(
  CMR* cmr,             /**< \brief \ref CMR environment. */
    CMR_SEYMOUR_NODE* dec  /**< \brief Decomposition node whose dense matrix to print. */
)
{
  CMR_UNUSED(cmr);
  CMR_UNUSED(dec);

#if defined(CMR_DEBUG)

  assert(dec);

  size_t* rowToMinorIndex = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowToMinorIndex, dec->numRows) );
  for (size_t row = 0; row < dec->numRows; ++row)
    rowToMinorIndex[row] = SIZE_MAX;

  size_t* columnToMinorIndex = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnToMinorIndex, dec->numColumns) );
  for (size_t column = 0; column < dec->numColumns; ++column)
    columnToMinorIndex[column] = SIZE_MAX;

  for (size_t m = 0; m < dec->nestedMinorsLength; ++m)
  {
    for (size_t r = (m == 0 ? 0 : dec->nestedMinorsSequenceNumRows[m - 1]);
      r < dec->nestedMinorsSequenceNumRows[m]; ++r)
    {
      rowToMinorIndex[dec->nestedMinorsRowsDense[r]] = m;
    }
    for (size_t c = (m == 0 ? 0 : dec->nestedMinorsSequenceNumColumns[m - 1]);
      c < dec->nestedMinorsSequenceNumColumns[m]; ++c)
    {
      columnToMinorIndex[dec->nestedMinorsColumnsDense[c]] = m;
    }
  }

  const char* letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  for (size_t row = 0; row < dec->numRows; ++row)
  {
    size_t rowMinor = rowToMinorIndex[row];
    for (size_t column = 0; column < dec->numColumns; ++column)
    {
      bool x = CMRdensebinmatrixGet(dec->denseMatrix, row, column);
      if (x)
      {
        size_t columnMinor = columnToMinorIndex[column];
        size_t minor = rowMinor > columnMinor ? rowMinor : columnMinor;
        putchar(minor == SIZE_MAX ? '1' : letters[minor % 26]);
      }
      else
        putchar(' ');
    }
    putchar('\n');
  }
  fflush(stdout);

  CMR_CALL( CMRfreeStackArray(cmr, &columnToMinorIndex) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowToMinorIndex) );

#endif /* CMR_DEBUG */
  return CMR_OKAY;
}

static
CMR_ERROR createHashVector(
  CMR* cmr,                 /**< \ref CMR environment. */
  long long** phashVector,  /**< Pointer for storing the hash vector. */
  size_t size               /**< Size of hash vector. */
)
{
  assert(cmr);

  CMR_CALL( CMRallocStackArray(cmr, phashVector, size) );
  long long* hashVector = *phashVector;
  size_t h = 1;
  for (size_t e = 0; e < size; ++e)
  {
    hashVector[e] = h;
    h = projectSignedHash(3 * h);
  }

  return CMR_OKAY;
}

/**
 * \brief Initializes the hashtable \p majorHashtable for unprocessed rows/columns.
 *
 * Initializes the hashtable \p majorHashtable for row/column vectors of unprocessed rows/columns with respect to those
 * columns/rows that were already processed.
 */

static
CMR_ERROR initializeHashing(
  CMR* cmr,                           /**< \ref CMR environment. */
  DenseBinaryMatrix* dense,           /**< Matrix. */
  ElementData* majorData,             /**< Major index data. */
  CMR_LISTHASHTABLE* majorHashtable,  /**< Major index hashtable. */
  size_t numMajors,                   /**< Number of major indices. */
  size_t* processedMinors,            /**< Array with processed minor indices. */
  size_t numProcessedMinors,          /**< Number of processed minor indices. */
  long long* hashVector,              /**< Hash vector. */
  bool isRow                          /**< Whether major means row. */
)
{
  assert(cmr);
  assert(dense);
  assert(majorData);
  assert(majorHashtable);
  assert(processedMinors);
  assert(hashVector);

  CMRdbgMsg(6, "Initializing %s hashtable.\n", isRow ? "row" : "column");

  /* Compute hash values for majors. */
  for (size_t major = 0; major < numMajors; ++major)
  {
    for (size_t j = 0; j < numProcessedMinors; ++j)
    {
      size_t minor = processedMinors[j];
      if (CMRdensebinmatrixGet(dense, isRow ? major : minor, isRow ? minor : major))
      {
        majorData[major].hashValue = projectSignedHash( majorData[major].hashValue + hashVector[minor] );
        majorData[major].numNonzeros++;
        majorData[major].representative = isRow ? CMRcolumnToElement(minor) : CMRrowToElement(minor);
      }
    }
    if (majorData[major].isProcessed)
    {
      CMR_CALL( CMRlisthashtableInsert(cmr, majorHashtable, majorData[major].hashValue, major,
        &majorData[major].hashEntry) );
    }
    else
      majorData[major].hashEntry = SIZE_MAX;
  }

  return CMR_OKAY;
}

/**
 * \brief Ensures that the representative entry for row/column \p index is correct.
 * 
 * Checks the respective hashtable.
 */

static
CMR_ERROR updateRepresentative(
  CMR* cmr,                           /**< \ref CMR environment. */
  DenseBinaryMatrix* dense,           /**< Matrix. */
  ElementData* majorData,             /**< Major index data. */
  CMR_LISTHASHTABLE* majorHashtable,  /**< Major index hashtable. */
  size_t* processedMinors,            /**< Array of processed minor indices. */
  size_t numProcessedMinors,          /**< Number of processed minor indices. */
  size_t majorIndex,                  /**< Major index to update. */
  bool isRow                          /**< Whether major means rows. */
)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(dense);
  assert(majorData);
  assert(majorHashtable);
  assert(processedMinors);

  if (majorData[majorIndex].isProcessed || majorData[majorIndex].numNonzeros <= 1)
    return CMR_OKAY;

  majorData[majorIndex].representative = 0;
  for (CMR_LISTHASHTABLE_ENTRY entry = CMRlisthashtableFindFirst(majorHashtable, majorData[majorIndex].hashValue);
    entry != SIZE_MAX; entry = CMRlisthashtableFindNext(majorHashtable, majorData[majorIndex].hashValue, entry))
  {
    size_t representativeIndex = CMRlisthashtableValue(majorHashtable, entry);
    if (majorData[majorIndex].hashValue != majorData[representativeIndex].hashValue)
      continue;

    bool equal = true;
    for (size_t j = 0; j < numProcessedMinors; ++j)
    {
      size_t index2 = processedMinors[j];
      bool e = CMRdensebinmatrixGet(dense, isRow ? majorIndex : index2, isRow ? index2 : majorIndex);
      bool f = CMRdensebinmatrixGet(dense, isRow ? representativeIndex : index2, isRow ? index2 : representativeIndex);
      if (e != f)
      {
        equal = false;
        break;
      }
    }
    if (equal)
    {
      majorData[majorIndex].representative = isRow ? CMRrowToElement(representativeIndex)
        : CMRcolumnToElement(representativeIndex);
      break;
    }
  }      

  return CMR_OKAY;
}

static
CMR_ERROR updateHashtable(
  CMR* cmr,                         /**< \ref CMR environment. */
  ElementData* majorData,           /**< Major index data. */
  size_t* processedMajors,          /**< Array of processed major indices. */
  size_t numProcessedMajors,        /**< Number of processed major indices. */
  CMR_LISTHASHTABLE* majorHashtable /**< Major index hashtable. */
)
{
  assert(cmr);
  assert(majorData);
  assert(processedMajors);
  assert(majorHashtable);

  CMRdbgMsg(8, "Updating row/column hashtable.\n");

  /* Add missing hashtable entries. */
  for (size_t i = 0; i < numProcessedMajors; i++)
  {
    size_t major = processedMajors[i];
    if (majorData[major].hashEntry == SIZE_MAX)
    {
      CMR_CALL( CMRlisthashtableInsert(cmr, majorHashtable, majorData[major].hashValue, major,
        &majorData[major].hashEntry) );
    }
  }

  return CMR_OKAY;
}

static
CMR_ERROR prepareSearch(
  CMR* cmr,                   /**< \ref CMR environment. */
  DenseBinaryMatrix* matrix,  /**< Nested matrix. */
  ElementData* rowData,       /**< Row data. */
  ElementData* columnData,    /**< Column data. */
  CMR_ELEMENT* pstartElement  /**< Pointer for storing an element for starting the search. */
)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(matrix);
  assert(rowData);
  assert(columnData);
  assert(pstartElement);

  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;

  CMRdbgMsg(6, "Preparing search.\n");

  for (size_t row = 0; row < numRows; ++row)
  {
    rowData[row].isSource = false;
    rowData[row].isTarget = false;
    rowData[row].isFlipped = false;
    if (!rowData[row].isProcessed && CMRelementIsValid(rowData[row].representative))
    {
      if (!*pstartElement)
        *pstartElement = rowData[row].representative;

      if (rowData[row].representative == *pstartElement)
      {
        rowData[row].isSource = true;
        if (CMRelementIsRow(*pstartElement))
          rowData[row].isFlipped = true;
      }
      else
        rowData[row].isTarget = true;
    }
  }
  for (size_t column = 0; column < numColumns; ++column)
  {
    columnData[column].isSource = false;
    columnData[column].isTarget = false;
    columnData[column].isFlipped = false;
    if (!columnData[column].isProcessed && CMRelementIsValid(columnData[column].representative))
    {
      if (!*pstartElement)
        *pstartElement = columnData[column].representative;
      if (columnData[column].representative == *pstartElement)
      {
        columnData[column].isSource = true;
        if (CMRelementIsColumn(*pstartElement))
          columnData[column].isFlipped = true;
      }
      else
        columnData[column].isTarget = true;
    }
  }

  if (CMRelementIsRow(*pstartElement))
  {
    size_t row = CMRelementToRowIndex(*pstartElement);
    for (size_t column = 0; column < matrix->numColumns; ++column)
    {
      if (CMRdensebinmatrixGet(matrix, row, column) && columnData[column].isTarget)
        columnData[column].isFlipped = true;
    }
  }
  else if (CMRelementIsColumn(*pstartElement))
  {
    size_t column = CMRelementToColumnIndex(*pstartElement);
    for (size_t row = 0; row < matrix->numRows; ++row)
    {
      if (CMRdensebinmatrixGet(matrix, row, column) && rowData[row].isTarget)
        rowData[row].isFlipped = true;
    }
  }

  return CMR_OKAY;
}


static
CMR_ERROR searchShortestPath(
  CMR* cmr,                   /**< \ref CMR environment. */
  DenseBinaryMatrix* dense,   /**< Remaining matrix. */
  ElementData* rowData,       /**< Row data. */
  ElementData* columnData,    /**< Column data. */
  CMR_ELEMENT* preachedTarget /**< Pointer for storing the reached target row/column. */
)
{
  assert(cmr);
  assert(dense);
  assert(rowData);
  assert(columnData);
  assert(preachedTarget);

  CMRdbgMsg(6, "Searching for shortest path.\n");
  *preachedTarget = 0;

  CMR_ELEMENT* queue = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &queue, dense->numRows + dense->numColumns) );
  size_t queueFirst = 0;
  size_t queueBeyond = 0;
  for (size_t row = 0; row < dense->numRows; ++row)
  {
    rowData[row].inQueue = rowData[row].isSource;
    if (rowData[row].isSource)
      queue[queueBeyond++] = CMRrowToElement(row);
    rowData[row].predecessor = 0;
  }
  for (size_t column = 0; column < dense->numColumns; ++column)
  {
    columnData[column].inQueue = columnData[column].isSource;
    if (columnData[column].isSource)
      queue[queueBeyond++] = CMRcolumnToElement(column);
    columnData[column].predecessor = 0;
  }

  while (queueFirst < queueBeyond)
  {
    CMRdbgMsg(8, "BFS queue contains %ld elements. Top element is %s.\n", queueBeyond - queueFirst,
    CMRelementString(queue[queueFirst], NULL));

    CMR_ELEMENT top = queue[queueFirst++];
    if (CMRelementIsRow(top))
    {
      size_t row = CMRelementToRowIndex(top);
      if (rowData[row].isTarget)
      {
        *preachedTarget = CMRrowToElement(row);
        break;
      }
      else
      {
        bool rowFlipped = rowData[row].isFlipped;
        for (size_t column = 0; column < dense->numColumns; ++column)
        {
          bool flip = rowFlipped && columnData[column].isFlipped;
          if (flip == CMRdensebinmatrixGet(dense, row, column))
            continue;
          if (columnData[column].inQueue || columnData[column].isProcessed)
            continue;
  
          queue[queueBeyond++] = CMRcolumnToElement(column);
          columnData[column].inQueue = true;
          columnData[column].predecessor = CMRrowToElement(row);
        }
      }
    }
    else
    {
      size_t column = CMRelementToColumnIndex(top);
      if (columnData[column].isTarget)
      {
        *preachedTarget = CMRcolumnToElement(column);
        break;
      }
      else
      {
        bool columnFlipped = columnData[column].isFlipped;
        for (size_t row = 0; row < dense->numRows; ++row)
        {
          bool flip = columnFlipped && rowData[row].isFlipped;
          if (flip == CMRdensebinmatrixGet(dense, row, column))
            continue;
          if (rowData[row].inQueue || rowData[row].isProcessed)
            continue;
  
          queue[queueBeyond++] = CMRrowToElement(row);
          rowData[row].inQueue = true;
          rowData[row].predecessor = CMRcolumnToElement(column);
        }
      }
    }
  }

  CMR_CALL( CMRfreeStackArray(cmr, &queue) );

  return CMR_OKAY;
}

/**
 * \brief Carries out a pivot in the dense matrix of \p dec.
 *
 * Also swaps the element entries of \p dec->denseRowsOriginal and \p dec->denseColumnsOriginal.
 */

static
CMR_ERROR pivot(
  CMR* cmr,                 /**< \ref CMR environment. */
    CMR_SEYMOUR_NODE* dec,     /**< decomposition node. */
  size_t pivotRow,          /**< Pivot row. */
  size_t pivotColumn        /**< Pivot column. */
)
{
  assert(cmr);
  assert(dec);

  /* Collect rows to be modified. */
  size_t* rows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rows, dec->numRows) );
  size_t numRows = 0;
  for (size_t row = 0; row < dec->numRows; ++row)
  {
    if (row != pivotRow && CMRdensebinmatrixGet(dec->denseMatrix, row, pivotColumn))
      rows[numRows++] = row;
  }

  /* Collect rows to be modified. */
  size_t* columns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columns, dec->numColumns) );
  size_t numColumns = 0;
  for (size_t column = 0; column < dec->numColumns; ++column)
  {
    if (column != pivotColumn && CMRdensebinmatrixGet(dec->denseMatrix, pivotRow, column))
      columns[numColumns++] = column;
  }

  for (size_t r = 0; r < numRows; ++r)
  {
    size_t row = rows[r];
    for (size_t c = 0; c < numColumns; ++c)
    {
      size_t column = columns[c];
      CMRdensebinmatrixFlip(dec->denseMatrix, row, column);
    }
  }
  
  CMR_CALL( CMRfreeStackArray(cmr, &columns) );
  CMR_CALL( CMRfreeStackArray(cmr, &rows) );

  CMR_ELEMENT temp = dec->denseRowsOriginal[pivotRow];
  dec->denseRowsOriginal[pivotRow] = dec->denseColumnsOriginal[pivotColumn];
  dec->denseColumnsOriginal[pivotColumn] = temp;
  
  return CMR_OKAY;
}

static
CMR_ERROR applyPivots(
  CMR* cmr,                   /**< \ref CMR environment. */
    CMR_SEYMOUR_NODE* dec,       /**< Decomposition node. */
  ElementData* rowData,       /**< Row data. */
  ElementData* columnData,    /**< Column data. */
  CMR_ELEMENT reachedTarget,  /**< Reached target row/column. */
  CMR_ELEMENT* newElements    /**< Array of length 3 for storing the new elements (followed by 0s). */
)
{
  assert(cmr);
  assert(dec);
  assert(rowData);
  assert(columnData);
  assert(CMRelementIsValid(reachedTarget));

  newElements[0] = 0;
  newElements[1] = 0;
  newElements[2] = 0;
  newElements[3] = 0;
  size_t countInvolvedElements = 0;
  for (CMR_ELEMENT e = reachedTarget; CMRelementIsValid(e); )
  {
    if (countInvolvedElements > 0)
    {
      newElements[2] = newElements[1];
      newElements[1] = newElements[0];
      newElements[0] = e;

      if (CMRelementIsValid(newElements[2]))
      {
        CMRdbgMsg(8, "Pivot at %s", CMRelementString(newElements[2], 0));
        CMRdbgMsg(0, ",%s.\n", CMRelementString(newElements[1], 0));

        unsigned char rowElement = CMRelementIsRow(newElements[2]) ? 2 : 1;
        size_t pivotRow = CMRelementToRowIndex(newElements[rowElement]);
        size_t pivotColumn = CMRelementToColumnIndex(newElements[3-rowElement]);
        CMR_CALL( pivot(cmr, dec, pivotRow, pivotColumn) );

        newElements[2] = 0;
        newElements[1] = 0;
      }
    }
    if (CMRelementIsRow(e))
      e = rowData[CMRelementToRowIndex(e)].predecessor;
    else
      e = columnData[CMRelementToColumnIndex(e)].predecessor;
    ++countInvolvedElements;
  }
  if (CMRelementIsValid(newElements[1]))
    newElements[2] = reachedTarget;
  else
    newElements[1] = reachedTarget;

  CMRdbgMsg(8, "Path involves %ld elements. Applied %ld pivots.\n", countInvolvedElements,
    (countInvolvedElements - 2) / 2);

  if (CMRelementIsValid(newElements[0]))
    CMRdbgMsg(10, "New element %s\n", CMRelementString(newElements[0], 0));
  if (CMRelementIsValid(newElements[1]))
    CMRdbgMsg(10, "New element %s\n", CMRelementString(newElements[1], 0));
  if (CMRelementIsValid(newElements[2]))
    CMRdbgMsg(10, "New element %s\n", CMRelementString(newElements[2], 0));

  return CMR_OKAY;
}

static
CMR_ERROR addElement(
  CMR* cmr,                           /**< \ref CMR environment. */
    CMR_SEYMOUR_NODE* dec,               /**< Decomposition node. */
  ElementData* majorData,             /**< Major index data. */
  ElementData* minorData,             /**< Minor index data. */
  CMR_LISTHASHTABLE* minorHashtable,  /**< Minor index hashtable. */
  long long* hashVector,              /**< Hash vector. */
  size_t numMinor,                    /**< Number of minor indices. */
  size_t* processedMajors,            /**< Array of processed major indices. */
  size_t* pnumProcessedMajors,        /**< Pointer to number of processed major indices. */
  size_t newMajor,                    /**< Major index to add. */
  size_t* nestedMinorsElements,       /**< Mapping from nested minor elements to elements of matrix. */
  bool isRow                          /**< Whether major means row. */
)
{
  assert(cmr);
  assert(dec);
  assert(majorData);
  assert(minorData);
  assert(minorHashtable);
  assert(!majorData[newMajor].isProcessed);

  CMRdbgMsg(8, "Adding %c%zu to processed submatrix.\n", isRow ? 'r' : 'c', newMajor+1);

  majorData[newMajor].isProcessed = true;
  processedMajors[*pnumProcessedMajors] = newMajor;
  (*pnumProcessedMajors)++;

  for (size_t minor = 0; minor < numMinor; ++minor)
  {
    if (CMRdensebinmatrixGet(dec->denseMatrix, isRow ? newMajor : minor, isRow ? minor : newMajor))
    {
      if (minorData[minor].isProcessed)
      {
        /* Remove minor from the hashtable (re-insertion after adding all elements). */
        if (minorData[minor].hashEntry != SIZE_MAX)
        {
          CMR_CALL( CMRlisthashtableRemove(cmr, minorHashtable, minorData[minor].hashEntry) );
          minorData[minor].hashEntry = SIZE_MAX;
        }
      }
      else
      {
        /* Increment nonzero counter. Store major as representative to take care of unit vectors. */
        minorData[minor].numNonzeros++;
        minorData[minor].representative = isRow ? CMRrowToElement(newMajor) : CMRcolumnToElement(newMajor);
      }
      /* In any case, update the hash value. */
      minorData[minor].hashValue = projectSignedHash( minorData[minor].hashValue + hashVector[newMajor] );
    }
  }

  /* Add to nestedMinorsRows or nestedMinorsColumns. */
  if (isRow)
    nestedMinorsElements[dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength-1]++] = newMajor;
  else
    nestedMinorsElements[dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength-1]++] = newMajor;

  return CMR_OKAY;
}

CMR_ERROR CMRregularityExtendNestedMinorSequence(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(queue);

    CMR_SEYMOUR_NODE* dec = task->node;
  assert(dec);

  CMRdbgMsg(6, "Attempting to extend a sequence of 3-connected nested minors of length %zu with "
    "last minor of size %zux%zu with the following dense matrix:\n",
    dec->nestedMinorsLength, dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength-1],
    dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength-1]);
  CMR_CALL( dbgPrintDenseSequence(cmr, dec) );

  size_t numRows = dec->matrix->numRows;
  size_t numColumns = dec->matrix->numColumns;
  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector, numRows > numColumns ? numRows : numColumns) );

  /* Initialize row data. */
  ElementData* rowData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowData, numRows) );
  for (size_t row = 0; row < numRows; ++row)
  {
    rowData[row].hashValue = 0;
    rowData[row].hashEntry = SIZE_MAX;
    rowData[row].representative = 0;
    rowData[row].numNonzeros = 0;
    rowData[row].isProcessed = false;
  }

  /* Initialize column data. */
  ElementData* columnData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnData, dec->matrix->numColumns) );
  for (size_t column = 0; column < numColumns; ++column)
  {
    columnData[column].hashValue = 0;
    columnData[column].hashEntry = SIZE_MAX;
    columnData[column].representative = 0;
    columnData[column].numNonzeros = 0;
    columnData[column].isProcessed = false;
  }

  /* Initialize information about existing sequence. */
  size_t* processedRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &processedRows, numRows) );
  size_t numProcessedRows = dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength-1];
  for (size_t r = 0; r < numProcessedRows; ++r)
  {
    size_t row = dec->nestedMinorsRowsDense[r];
    rowData[row].representative = CMRrowToElement(row);
    rowData[row].isProcessed = true;
    processedRows[r] = row;
  }

  size_t* processedColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &processedColumns, numColumns) );
  size_t numProcessedColumns = dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength-1];
  for (size_t c = 0; c < numProcessedColumns; ++c)
  {
    size_t column = dec->nestedMinorsColumnsDense[c];
    columnData[column].representative = CMRcolumnToElement(column);
    columnData[column].isProcessed = true;
    processedColumns[c] = column;
  }

  /* Create hash maps. */
  CMR_LISTHASHTABLE* rowHashtable = NULL;
  CMR_CALL( CMRlisthashtableCreate(cmr, &rowHashtable, nextPower2(numRows), numRows) );
  CMR_LISTHASHTABLE* columnHashtable = NULL;
  CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns) );

  CMR_CALL( initializeHashing(cmr, dec->denseMatrix, rowData, rowHashtable, numRows, processedColumns,
    numProcessedColumns, hashVector, true) );
  CMR_CALL( initializeHashing(cmr, dec->denseMatrix, columnData, columnHashtable, numColumns, processedRows,
    numProcessedRows, hashVector, false) );

  CMR_ERROR result = CMR_OKAY;
  size_t elementTimeFactor = (numRows + numColumns) / 100 + 1;
  while (numProcessedRows < numRows || numProcessedColumns < numColumns)
  {
    if (((numProcessedRows + numProcessedColumns) % elementTimeFactor == 0)
      && (clock() - task->startClock) * 1.0 / CLOCKS_PER_SEC > task->timeLimit)
    {
      goto cleanup;
    }

    CMRdbgMsg(8, "New iteration; processed %zu rows and %zu columns so far.\n", numProcessedRows, numProcessedColumns);
    if (task->stats)
      task->stats->sequenceExtensionCount++;

    CMR_CALL( dbgPrintDenseSequence(cmr, dec) );

    for (size_t row = 0; row < numRows; ++row)
    {
      if (rowData[row].isProcessed)
        CMRdbgMsg(10, "Row r%zu was already processed.\n", row+1);
      else if (rowData[row].numNonzeros == 0)
        CMRdbgMsg(10, "Row r%zu is a zero row.\n", row+1);
      else if (rowData[row].numNonzeros == 1)
        CMRdbgMsg(10, "Row r%zu is a unit row for element %s.\n", row+1,
          CMRelementString(rowData[row].representative, 0));
      else
        CMRdbgMsg(10, "Row r%zu may be parallel.\n", row+1);
    }

    for (size_t column = 0; column < numColumns; ++column)
    {
      if (columnData[column].isProcessed)
        CMRdbgMsg(10, "Column c%zu was already processed.\n", column+1);
      else if (columnData[column].numNonzeros == 0)
        CMRdbgMsg(10, "Column c%zu is a zero column.\n", column+1);
      else if (columnData[column].numNonzeros == 1)
        CMRdbgMsg(10, "Column c%zu is a unit row for element %s.\n", column+1,
          CMRelementString(columnData[column].representative, 0));
      else
        CMRdbgMsg(10, "Column c%zu may be parallel.\n", column+1);
    }

    bool added = false;
    for (size_t row = 0; row < numRows; ++row)
    {
      if (!rowData[row].isProcessed && rowData[row].numNonzeros > 1)
      {
        CMR_CALL( updateRepresentative(cmr, dec->denseMatrix, rowData, rowHashtable, processedColumns,
          numProcessedColumns, row, true) );
        if (CMRelementIsValid(rowData[row].representative))
        {
          CMRdbgMsg(10, "Row r%zu is parallel to processed row %s\n", row+1,
            CMRelementString(rowData[row].representative, 0));
        }
        else
        {
          CMRdbgMsg(10, "Encountered non-parallel row r%zu.\n", row+1);
          dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength]
            = dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength-1];
          dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength]
            = dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength-1];
          dec->nestedMinorsLength++;
          CMR_CALL( addElement(cmr, dec, rowData, columnData, columnHashtable, hashVector, numColumns, processedRows,
            &numProcessedRows, row, dec->nestedMinorsRowsDense, true) );
          CMR_CALL( updateHashtable(cmr, rowData, processedRows, numProcessedRows, rowHashtable) );
          CMR_CALL( updateHashtable(cmr, columnData, processedColumns, numProcessedColumns, columnHashtable) );
          added = true;
          break;
        }
      }
    }
    if (added)
      continue;

    for (size_t column = 0; column < numColumns; ++column)
    {
      if (!columnData[column].isProcessed && columnData[column].numNonzeros > 1)
      {
        CMR_CALL( updateRepresentative(cmr, dec->denseMatrix, columnData, columnHashtable, processedRows,
          numProcessedRows, column, false) );
        if (CMRelementIsValid(columnData[column].representative))
        {
          CMRdbgMsg(8, "Column c%zu is parallel to processed column %s\n", column+1,
            CMRelementString(columnData[column].representative, 0));
        }
        else
        {
          CMRdbgMsg(8, "Encountered non-parallel column c%zu.\n", column+1);
          dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength]
            = dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength-1];
          dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength]
          = dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength-1];
          dec->nestedMinorsLength++;
          CMR_CALL( addElement(cmr, dec, columnData, rowData, rowHashtable, hashVector, numRows, processedColumns,
            &numProcessedColumns, column, dec->nestedMinorsColumnsDense, false) );
          CMR_CALL( updateHashtable(cmr, rowData, processedRows, numProcessedRows, rowHashtable) );
          CMR_CALL( updateHashtable(cmr, columnData, processedColumns, numProcessedColumns, columnHashtable) );
          added = true;
          break;
        }
      }
    }
    if (added)
      continue;

    /* All unprocessed rows/columns are zero, unit or parallel to submatrix. */
    CMR_ELEMENT startElement = 0;
    CMR_CALL( prepareSearch(cmr, dec->denseMatrix, rowData, columnData, &startElement) );
    for (size_t row = 0; row < numRows; ++row)
    {
      if (rowData[row].isSource)
        CMRdbgMsg(10, "Row r%zu is a source.\n", row+1);
      if (rowData[row].isTarget)
        CMRdbgMsg(10, "Row r%zu is a target.\n", row+1);
      if (rowData[row].isFlipped)
        CMRdbgMsg(10, "Row r%zu is flipped.\n", row+1);
    }

    for (size_t column = 0; column < numColumns; ++column)
    {
      if (columnData[column].isSource)
        CMRdbgMsg(10, "Column c%zu is a source.\n", column+1);
      else if (columnData[column].isTarget)
        CMRdbgMsg(10, "Column c%zu is a target.\n", column+1);
      else if (columnData[column].isFlipped)
        CMRdbgMsg(10, "Column c%zu is flipped.\n", column+1);
    }

    CMR_ELEMENT reachedTarget;
    CMR_CALL( searchShortestPath(cmr, dec->denseMatrix, rowData, columnData, &reachedTarget) );

    if (CMRelementIsValid(reachedTarget))
    {
      CMRdbgMsg(10, "Path exists:\n");
      for (CMR_ELEMENT e = reachedTarget; CMRelementIsValid(e); )
      {
        CMRdbgMsg(12, "%s\n", CMRelementString(e, 0));
        if (CMRelementIsRow(e))
          e = rowData[CMRelementToRowIndex(e)].predecessor;
        else
          e = columnData[CMRelementToColumnIndex(e)].predecessor;
      }

      dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength]
        = dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength-1];
      dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength]
        = dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength-1];
      dec->nestedMinorsLength++;
      CMR_ELEMENT newElements[4];
      CMR_CALL( applyPivots(cmr, dec, rowData, columnData, reachedTarget, newElements) );

      for (size_t i = 0; CMRelementIsValid(newElements[i]); ++i)
      {
        if (CMRelementIsRow(newElements[i]))
        {
          CMR_CALL( addElement(cmr, dec, rowData, columnData, columnHashtable, hashVector, numColumns, processedRows,
            &numProcessedRows, CMRelementToRowIndex(newElements[i]), dec->nestedMinorsRowsDense, true) );
        }
        else
        {
          CMR_CALL( addElement(cmr, dec, columnData, rowData, rowHashtable, hashVector, numRows, processedColumns,
            &numProcessedColumns, CMRelementToColumnIndex(newElements[i]), dec->nestedMinorsColumnsDense, false) );
        }
      }

      CMR_CALL( updateHashtable(cmr, rowData, processedRows, numProcessedRows, rowHashtable) );
      CMR_CALL( updateHashtable(cmr, columnData, processedColumns, numProcessedColumns, columnHashtable) );
    }
    else
    {
      CMRdbgMsg(10, "No path from start element %s found.\n", CMRelementString(startElement, 0));

      /* Create the separation object. */
      CMR_SEPA* separation = NULL;
      CMR_CALL( CMRsepaCreate(cmr, numRows, numColumns, &separation) );

      for (size_t row = 0; row < numRows; ++row)
      {
        CMR_SEPA_FLAGS flag = (CMRelementIsValid(rowData[row].predecessor) || rowData[row].isSource
          || CMRrowToElement(row) == startElement) ? CMR_SEPA_SECOND : CMR_SEPA_FIRST;
        CMRdbgMsg(12, "Row r%zu has predecessor %s, isSource=%s and lies in part %d.\n", row+1,
          CMRelementString(rowData[row].predecessor, 0), rowData[row].isSource ? "true" : "false", flag);
        CMR_ELEMENT originalElement = dec->denseRowsOriginal[row];
        if (CMRelementIsRow(originalElement))
          separation->rowsFlags[CMRelementToRowIndex(originalElement)] = flag;
        else
          separation->columnsFlags[CMRelementToColumnIndex(originalElement)] = flag;
      }

      for (size_t column = 0; column < numColumns; ++column)
      {
        CMR_SEPA_FLAGS flag = (CMRelementIsValid(columnData[column].predecessor) || columnData[column].isSource
          || CMRcolumnToElement(column) == startElement) ? CMR_SEPA_SECOND : CMR_SEPA_FIRST;
        CMRdbgMsg(12, "Column c%zu has predecessor %s, isSource=%s and lies in part %d.\n", column+1,
          CMRelementString(columnData[column].predecessor, 0), columnData[column].isSource ? "true" : "false", flag);
        CMR_ELEMENT originalElement = dec->denseColumnsOriginal[column];
        if (CMRelementIsRow(originalElement))
          separation->rowsFlags[CMRelementToRowIndex(originalElement)] = flag;
        else
          separation->columnsFlags[CMRelementToColumnIndex(originalElement)] = flag;
      }

      if (dec->transpose == NULL)
        CMR_CALL( CMRchrmatTranspose(cmr, dec->matrix, &dec->transpose) );

      CMR_SUBMAT* violatorSubmatrix = NULL;
      CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, separation, dec->matrix, dec->transpose, NULL,
        dec->isTernary ? &violatorSubmatrix : NULL) );

      if (violatorSubmatrix)
      {
        CMRdbgMsg(8, "-> 2x2 submatrix with bad determinant.\n");

        CMR_CALL( CMRseymourUpdateViolator(cmr, dec, violatorSubmatrix) );
        assert(dec->type == CMR_SEYMOUR_NODE_TYPE_IRREGULAR);

        CMR_CALL( CMRsepaFree(cmr, &separation) );
        queue->foundIrregularity = true;

        break;
      }

      assert(separation->type == CMR_SEPA_TYPE_TWO);

      /* Carry out 2-sum decomposition. */

      CMRdbgMsg(8, "-> 2-separation found.\n");
      CMR_CALL( CMRseymourUpdateTwoSum(cmr, dec, separation) );
      CMR_CALL( CMRsepaFree(cmr, &separation) );

      DecompositionTask* childTasks[2] = { task, NULL };
      CMR_CALL( CMRregularityTaskCreateRoot(cmr, dec->children[1], &childTasks[1], task->params, task->stats,
        task->startClock, task->timeLimit) );

      childTasks[0]->node = dec->children[0];
      dec->children[0]->testedSeriesParallel = false; /* TODO: we may carry over the found sequence including W_k. */
      dec->children[1]->testedSeriesParallel = false;

      /* Add both child tasks to the list. */
      CMRregularityQueueAdd(queue, childTasks[0]);
      CMRregularityQueueAdd(queue, childTasks[1]);

      break;
    }
  }

  if (task->stats)
  {
    task->stats->sequenceExtensionTime += (clock() - task->startClock) * 1.0 / CLOCKS_PER_SEC;
  }

  if (dec->type == CMR_SEYMOUR_NODE_TYPE_TWO_SUM)
  {
    CMRdbgMsg(8, "Aborting construction of sequence of nested 3-connected minors due to a 2-separation.\n");
  }
  else if (dec->type == CMR_SEYMOUR_NODE_TYPE_IRREGULAR)
  {
    CMRdbgMsg(8, "Aborting construction of sequence of nested 3-connected minors due to detected irregularity.\n");

    /* Task is done. */
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }
  else
  {
    CMRdbgMsg(8, "Successfully constructed sequence of nested 3-connected minors.\n");
    assert(dec->type != CMR_SEYMOUR_NODE_TYPE_IRREGULAR);

    assert(dec->nestedMinorsMatrix == NULL);
    assert(dec->nestedMinorsTranspose == NULL);

    /* Count the nonzeros. */
    size_t entry = 0;
    for (size_t row = 0; row < numRows; ++row)
    {
      for (size_t column = 0; column < numColumns; ++column)
      {
        if (CMRdensebinmatrixGet(dec->denseMatrix, row, column))
          ++entry;
      }
    }

    /* Create a sparse copy of dense, permuted such that the nested sequence is displayed from top-left on. */
    CMR_CALL( CMRchrmatCreate(cmr, &dec->nestedMinorsMatrix, numRows, numColumns, entry) );
    entry = 0;
    for (size_t row = 0; row < numRows; ++row)
    {
      size_t denseRow = dec->nestedMinorsRowsDense[row];
      dec->nestedMinorsMatrix->rowSlice[row] = entry;
      for (size_t column = 0; column < numColumns; ++column)
      {
        size_t denseColumn = dec->nestedMinorsColumnsDense[column];
        if (CMRdensebinmatrixGet(dec->denseMatrix, denseRow, denseColumn))
        {
          dec->nestedMinorsMatrix->entryColumns[entry] = column;
          dec->nestedMinorsMatrix->entryValues[entry] = 1;
          ++entry;
        }
      }
    }
    dec->nestedMinorsMatrix->rowSlice[numRows] = entry;
    assert(entry == dec->nestedMinorsMatrix->numNonzeros);

    /* Create the transpose. */
    assert(dec->nestedMinorsTranspose == NULL);
    CMR_CALL( CMRchrmatTranspose(cmr, dec->nestedMinorsMatrix, &dec->nestedMinorsTranspose) );

    /* Create mappings for nested minors matrix. */
    assert(dec->nestedMinorsRowsOriginal == NULL);
    CMR_CALL( CMRallocBlockArray(cmr, &dec->nestedMinorsRowsOriginal, dec->numRows) );
    for (size_t row = 0; row < dec->numRows; ++row)
      dec->nestedMinorsRowsOriginal[row] = dec->denseRowsOriginal[dec->nestedMinorsRowsDense[row]];

    assert(dec->nestedMinorsColumnsOriginal == NULL);
    CMR_CALL( CMRallocBlockArray(cmr, &dec->nestedMinorsColumnsOriginal, dec->numColumns) );
    for (size_t column = 0; column < dec->numColumns; ++column)
      dec->nestedMinorsColumnsOriginal[column] = dec->denseColumnsOriginal[dec->nestedMinorsColumnsDense[column]];

    /* Free all the dense matrix data. */
    CMR_CALL( CMRdensebinmatrixFree(cmr, &dec->denseMatrix) );
    CMR_CALL( CMRfreeBlockArray(cmr, &dec->denseRowsOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &dec->denseColumnsOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsRowsDense) );
    CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsColumnsDense) );

    /* Add the task back to the list of unprocessed tasks in order to check how far the sequence is (co)graphic. */
    CMRregularityQueueAdd(queue, task);
  }

cleanup:

  CMR_CALL( CMRlisthashtableFree(cmr, &columnHashtable) );
  CMR_CALL( CMRlisthashtableFree(cmr, &rowHashtable) );

  CMR_CALL( CMRfreeStackArray(cmr, &processedColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &processedRows) );

  CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowData) );
  CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );

  return result;
}

CMR_ERROR CMRregularityInitNestedMinorSequence(CMR* cmr, DecompositionTask* task, CMR_SUBMAT* wheelSubmatrix)
{
  assert(cmr);
  assert(task);
  assert(wheelSubmatrix);

    CMR_SEYMOUR_NODE* dec = task->node;
  assert(dec);

#if defined(CMR_DEBUG)
  CMRdbgMsg(8, "Initializing a sequence of nested 3-connected minors for the following matrix:\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stdout, '0', true) );
  CMRdbgMsg(10, "Row indices: ");
  for (size_t row = 0; row < wheelSubmatrix->numRows; ++row)
    CMRdbgMsg(0, "%zu ", wheelSubmatrix->rows[row]);
  CMRdbgMsg(0, "\n");
  CMRdbgMsg(10, "\nColumn indices: ");
  for (size_t column = 0; column < wheelSubmatrix->numColumns; ++column)
    CMRdbgMsg(0, "%zu ", wheelSubmatrix->columns[column]);
  CMRdbgMsg(0, "\n");
#endif /* CMR_DEBUG */

  size_t numRows = dec->matrix->numRows;
  size_t numColumns = dec->matrix->numColumns;

  CMR_CALL( CMRdensebinmatrixCreate(cmr, numRows, numColumns, &dec->denseMatrix) );
  CMR_CALL( CMRallocBlockArray(cmr, &dec->nestedMinorsRowsDense, numRows) );
  CMR_CALL( CMRallocBlockArray(cmr, &dec->nestedMinorsColumnsDense, numColumns) );
  CMR_CALL( CMRallocBlockArray(cmr, &dec->denseRowsOriginal, numRows) );
  CMR_CALL( CMRallocBlockArray(cmr, &dec->denseColumnsOriginal, numColumns) );

  CMR_CALL( CMRallocBlockArray(cmr, &dec->nestedMinorsSequenceNumRows, numRows + numColumns) );
  CMR_CALL( CMRallocBlockArray(cmr, &dec->nestedMinorsSequenceNumColumns, numRows + numColumns) );

  /* Copy the matrix into the dense one. */
  for (size_t row = 0; row < numRows; ++row)
  {
    size_t first = dec->matrix->rowSlice[row];
    size_t beyond = dec->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
      CMRdensebinmatrixSet1(dec->denseMatrix, row, dec->matrix->entryColumns[e]);
  }

  /* Start keeping track of elements in order to trace back pivots. */
  for (size_t row = 0; row < numRows; ++row)
    dec->denseRowsOriginal[row] = CMRrowToElement(row);
  for (size_t column = 0; column < numColumns; ++column)
    dec->denseColumnsOriginal[column] = CMRcolumnToElement(column);

  dec->nestedMinorsLength = 1;
  dec->nestedMinorsSequenceNumRows[0] = 3;
  dec->nestedMinorsSequenceNumColumns[0] = 3;

  for (size_t i = 0; i < 3; ++i)
  {
    dec->nestedMinorsRowsDense[i] = wheelSubmatrix->rows[i];
    dec->nestedMinorsColumnsDense[i] = wheelSubmatrix->columns[i];
  }

  CMRdbgMsg(10, "Dense matrix before pivoting to W_3:\n");
  CMR_CALL( dbgPrintDenseSequence(cmr, dec) );

  /* Carry out pivots on dense matrix to shorten the W_k to W_3. */

  for (size_t p = wheelSubmatrix->numRows - 1; p >= 3; --p)
  {
    CMR_CALL( pivot(cmr, dec, wheelSubmatrix->rows[p], wheelSubmatrix->columns[p]) );

    CMRdbgMsg(10, "After pivoting at %zu,%zu.\n", wheelSubmatrix->rows[p], wheelSubmatrix->columns[p]);
    CMR_CALL( dbgPrintDenseSequence(cmr, dec) );
  }

  if (task->stats)
    task->stats->sequenceExtensionCount++;

  return CMR_OKAY;
}

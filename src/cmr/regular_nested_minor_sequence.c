// #define CMR_DEBUG /* Uncomment to debug this file. */

#include "dec_internal.h"
#include "regular_internal.h"
#include "env_internal.h"

#include "densematrix.h"
#include "hashtable.h"

typedef struct
{
  long long hashValue;              /**< \brief Hash value of this element. */
  size_t hashEntry;                 /**< \brief Entry in hashtable. */
  size_t numNonzeros;               /**< \brief Number of nonzeros in (part parallel to) processed submatrix. */
  CMR_ELEMENT representative;       /**< \brief Parallel element. */
  CMR_ELEMENT predecessor;          /**< \brief Predecessor row/column in BFS. */
  CMR_ELEMENT original;             /**< \brief Element of given matrix that this row/column represents. */
  bool isProcessed : 1;             /**< \brief Whether this row/column belongs to processed submatrix. */
  bool isSource : 1;                /**< \brief Whether this row/column is a source node in the BFS. */
  bool isTarget : 1;                /**< \brief Whether this row/column is a target node in the BFS. */
  bool isFlipped : 1;               /**< \brief Whether this row/column is a flip node. Entries in a flip row and a
                                                flip column change roles with respect to being edges. */
  bool inQueue : 1;                 /**< \brief Whether this row/column is in the BFS queue. */
} ElementData;

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
  DenseBinaryMatrix* dense,         /**< Matrix. */
  ElementData* majorData,           /**< Major index data. */
  size_t* processedMajors,          /**< Array of processed major indices. */
  size_t numProcessedMajors,        /**< Number of processed major indices. */
  CMR_LISTHASHTABLE* majorHashtable /**< Major index hashtable. */
)
{
  assert(cmr);
  assert(dense);
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

static
CMR_ERROR pivot(
  CMR* cmr,                 /**< \ref CMR environment. */
  DenseBinaryMatrix* dense, /**< Remaining matrix. */
  size_t pivotRow,          /**< Pivot row. */
  size_t pivotColumn        /**< Pivot column. */
)
{
  assert(cmr);
  assert(dense);

  /* Collect rows to be modified. */
  size_t* rows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rows, dense->numRows) );
  size_t numRows = 0;
  for (size_t row = 0; row < dense->numRows; ++row)
  {
    if (row != pivotRow && CMRdensebinmatrixGet(dense, row, pivotColumn))
      rows[numRows++] = row;
  }

  /* Collect rows to be modified. */
  size_t* columns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columns, dense->numColumns) );
  size_t numColumns = 0;
  for (size_t column = 0; column < dense->numColumns; ++column)
  {
    if (column != pivotColumn && CMRdensebinmatrixGet(dense, pivotRow, column))
      columns[numColumns++] = column;
  }

  for (size_t r = 0; r < numRows; ++r)
  {
    size_t row = rows[r];
    for (size_t c = 0; c < numColumns; ++c)
    {
      size_t column = columns[c];
      CMRdensebinmatrixFlip(dense, row, column);
    }
  }
  
  CMR_CALL( CMRfreeStackArray(cmr, &columns) );
  CMR_CALL( CMRfreeStackArray(cmr, &rows) );
  
  return CMR_OKAY;
}

static
CMR_ERROR applyPivots(
  CMR* cmr,                   /**< \ref CMR environment. */
  DenseBinaryMatrix* dense,   /**< Remaining matrix. */
  ElementData* rowData,       /**< Row data. */
  ElementData* columnData,    /**< Column data. */
  CMR_ELEMENT reachedTarget,  /**< Reached target row/column. */
  CMR_ELEMENT* newElements    /**< Array of length 3 for storing the new elements (followed by 0s). */
)
{
  assert(cmr);
  assert(dense);
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
        CMR_CALL( pivot(cmr, dense, pivotRow, pivotColumn) );

        /* Exchange links to original row/column. */
        CMR_ELEMENT tmp = rowData[pivotRow].original;
        rowData[pivotRow].original = columnData[pivotColumn].original;
        columnData[pivotColumn].original = tmp;

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
  DenseBinaryMatrix* dense,           /**< Matrix. */
  ElementData* majorData,             /**< Major index data. */
  ElementData* minorData,             /**< Minor index data. */
  CMR_LISTHASHTABLE* minorHashtable,  /**< Minor index hashtable. */
  long long* hashVector,              /**< Hash vector. */
  size_t numMinor,                    /**< Number of minor indices. */
  size_t* processedMajors,            /**< Array of processed major indices. */
  size_t* pnumProcessedMajors,        /**< Pointer to number of processed major indices. */
  size_t newMajor,                    /**< Major index to add. */
  bool isRow                          /**< Whether major means row. */
)
{
  assert(cmr);
  assert(dense);
  assert(majorData);
  assert(minorData);
  assert(minorHashtable);
  assert(!majorData[newMajor].isProcessed);

  CMRdbgMsg(8, "Adding %c%ld to processed submatrix.\n", isRow ? 'r' : 'c', newMajor+1);

  majorData[newMajor].isProcessed = true;
  processedMajors[*pnumProcessedMajors] = newMajor;
  (*pnumProcessedMajors)++;

  for (size_t minor = 0; minor < numMinor; ++minor)
  {
    if (CMRdensebinmatrixGet(dense, isRow ? newMajor : minor, isRow ? minor : newMajor))
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

  return CMR_OKAY;
}

CMR_ERROR CMRregularConstructNestedMinorSequence(CMR* cmr, CMR_DEC* dec, bool ternary, CMR_SUBMAT* wheelSubmatrix,
  CMR_SUBMAT** psubmatrix, CMR_REGULAR_PARAMETERS* params)
{
  assert(cmr);
  assert(dec);
  assert(wheelSubmatrix);
  assert(params);

  CMRdbgMsg(4, "Attempting to construct a sequence of 3-connected nested minors.\n");

  CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stdout, '0', true) );
  
  size_t numRows = dec->matrix->numRows;
  size_t numColumns = dec->matrix->numColumns;
  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector, numRows > numColumns ? numRows : numColumns) );

  /* Create dense matrix. */
  DenseBinaryMatrix* dense = NULL;
  CMR_CALL( CMRdensebinmatrixCreateStack(cmr, dec->matrix->numRows, dec->matrix->numColumns, &dense) );
  for (size_t row = 0; row < dense->numRows; ++row)
  {
    size_t first = dec->matrix->rowSlice[row];
    size_t beyond = dec->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
      CMRdensebinmatrixSet1(dense, row, dec->matrix->entryColumns[e]);
  }

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
    rowData[row].original = CMRrowToElement(row);
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
    columnData[column].original = CMRcolumnToElement(column);
  }

  /* Initialize already-processed part. */
  size_t* processedRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &processedRows, numRows) );
  size_t numProcessedRows = wheelSubmatrix->numRows;
  size_t* processedColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &processedColumns, numColumns) );
  size_t numProcessedColumns = wheelSubmatrix->numColumns;
  for (size_t r = 0; r < wheelSubmatrix->numRows; ++r)
  {
    size_t row = wheelSubmatrix->rows[r];
    rowData[row].representative = CMRrowToElement(row);
    rowData[row].isProcessed = true;
    processedRows[r] = row;
  }
  for (size_t c = 0; c < wheelSubmatrix->numColumns; ++c)
  {
    size_t column = wheelSubmatrix->columns[c];
    columnData[column].representative = CMRcolumnToElement(column);
    columnData[column].isProcessed = true;
    processedColumns[c] = column;
  }

  /* Create hash maps. */
  CMR_LISTHASHTABLE* rowHashtable = NULL;
  CMR_CALL( CMRlisthashtableCreate(cmr, &rowHashtable, nextPower2(numRows), numRows) );
  CMR_LISTHASHTABLE* columnHashtable = NULL;
  CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns) );

  CMR_CALL( initializeHashing(cmr, dense, rowData, rowHashtable, numRows, processedColumns, numProcessedColumns, 
    hashVector, true) );
  CMR_CALL( initializeHashing(cmr, dense, columnData, columnHashtable, numColumns, processedRows, numProcessedRows, 
    hashVector, false) );

  while (numProcessedRows < numRows || numProcessedColumns < numColumns)
  {
    CMRdbgMsg(6, "New iteration; processed %ld rows and %ld columns so far.\n", numProcessedRows, numProcessedColumns);
    
    CMRdbgMsg(6, "Dense matrix:\n");
    for (size_t row = 0; row < dense->numRows; ++row)
    {
      CMRdbgMsg(8, "");
      for (size_t column = 0; column < dense->numColumns; ++column)
        CMRdbgMsg(0, " %d", CMRdensebinmatrixGet(dense, row, column));
      CMRdbgMsg(0, "\n");
    }
    
    for (size_t row = 0; row < numRows; ++row)
    {
      if (rowData[row].isProcessed)
        CMRdbgMsg(8, "Row r%ld was already processed.\n", row+1);
      else if (rowData[row].numNonzeros == 0)
        CMRdbgMsg(8, "Row r%ld is a zero row.\n", row+1);
      else if (rowData[row].numNonzeros == 1)
        CMRdbgMsg(8, "Row r%ld is a unit row for element %s.\n", row+1,
          CMRelementString(rowData[row].representative, 0));
      else
        CMRdbgMsg(8, "Row r%ld may be parallel.\n", row+1);
    }

    for (size_t column = 0; column < numColumns; ++column)
    {
      if (columnData[column].isProcessed)
        CMRdbgMsg(8, "Column c%ld was already processed.\n", column+1);
      else if (columnData[column].numNonzeros == 0)
        CMRdbgMsg(8, "Column c%ld is a zero column.\n", column+1);
      else if (columnData[column].numNonzeros == 1)
        CMRdbgMsg(8, "Column c%ld is a unit row for element %s.\n", column+1,
          CMRelementString(columnData[column].representative, 0));
      else
        CMRdbgMsg(8, "Column c%ld may be parallel.\n", column+1);
    }

    bool added = false;
    for (size_t row = 0; row < numRows; ++row)
    {
      if (!rowData[row].isProcessed && rowData[row].numNonzeros > 1)
      {
        CMR_CALL( updateRepresentative(cmr, dense, rowData, rowHashtable, processedColumns, numProcessedColumns, row,
          true) );
        if (CMRelementIsValid(rowData[row].representative))
        {
          CMRdbgMsg(8, "Row r%ld is parallel to processed row %s\n", row+1,
            CMRelementString(rowData[row].representative, 0));
        }
        else
        {
          CMRdbgMsg(8, "Encountered non-parallel row r%ld.\n", row+1);
          CMR_CALL( addElement(cmr, dense, rowData, columnData, columnHashtable, hashVector, numColumns, processedRows,
            &numProcessedRows, row, true) );
          CMR_CALL( updateHashtable(cmr, dense, rowData, processedRows, numProcessedRows, rowHashtable) );
          CMR_CALL( updateHashtable(cmr, dense, columnData, processedColumns, numProcessedColumns, columnHashtable) );
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
        CMR_CALL( updateRepresentative(cmr, dense, columnData, columnHashtable, processedRows, numProcessedRows, column,
          false) );
        if (CMRelementIsValid(columnData[column].representative))
        {
          CMRdbgMsg(8, "Column c%ld is parallel to processed column %s\n", column+1,
            CMRelementString(columnData[column].representative, 0));
        }
        else
        {
          CMRdbgMsg(8, "Encountered non-parallel column c%ld.\n", column+1);
          CMR_CALL( addElement(cmr, dense, columnData, rowData, rowHashtable, hashVector, numRows, processedColumns,
            &numProcessedColumns, column, false) );
          CMR_CALL( updateHashtable(cmr, dense, rowData, processedRows, numProcessedRows, rowHashtable) );
          CMR_CALL( updateHashtable(cmr, dense, columnData, processedColumns, numProcessedColumns, columnHashtable) );
          added = true;
          break;
        }
      }
    }
    if (added)
      continue;

    /* All unprocessed rows/columns are zero, unit or parallel to submatrix. */

    CMR_ELEMENT startElement = 0;
    CMR_CALL( prepareSearch(cmr, dense, rowData, columnData, &startElement) );

    for (size_t row = 0; row < numRows; ++row)
    {
      if (rowData[row].isSource)
        CMRdbgMsg(8, "Row r%ld is a source.\n", row+1);
      if (rowData[row].isTarget)
        CMRdbgMsg(8, "Row r%ld is a target.\n", row+1);
      if (rowData[row].isFlipped)
        CMRdbgMsg(8, "Row r%ld is flipped.\n", row+1);
    }

    for (size_t column = 0; column < numColumns; ++column)
    {
      if (columnData[column].isSource)
        CMRdbgMsg(8, "Column c%ld is a source.\n", column+1);
      else if (columnData[column].isTarget)
        CMRdbgMsg(8, "Column c%ld is a target.\n", column+1);
      else if (columnData[column].isFlipped)
        CMRdbgMsg(8, "Column c%ld is flipped.\n", column+1);
    }

    CMR_ELEMENT reachedTarget;
    CMR_CALL( searchShortestPath(cmr, dense, rowData, columnData, &reachedTarget) );

    if (CMRelementIsValid(reachedTarget))
    {
      CMRdbgMsg(8, "Path exists:\n");
      for (CMR_ELEMENT e = reachedTarget; CMRelementIsValid(e); )
      {
        CMRdbgMsg(10, "%s\n", CMRelementString(e, 0));
        if (CMRelementIsRow(e))
          e = rowData[CMRelementToRowIndex(e)].predecessor;
        else
          e = columnData[CMRelementToColumnIndex(e)].predecessor;
      }

      CMR_ELEMENT newElements[4];
      CMR_CALL( applyPivots(cmr, dense, rowData, columnData, reachedTarget, newElements) );

      for (size_t i = 0; CMRelementIsValid(newElements[i]); ++i)
      {
        if (CMRelementIsRow(newElements[i]))
        {
          CMR_CALL( addElement(cmr, dense, rowData, columnData, columnHashtable, hashVector, numColumns,
            processedRows, &numProcessedRows, CMRelementToRowIndex(newElements[i]), true) );
        }
        else
        {
          CMR_CALL( addElement(cmr, dense, columnData, rowData, rowHashtable, hashVector, numRows,
            processedColumns, &numProcessedColumns, CMRelementToColumnIndex(newElements[i]), false) );
        }
      }

      CMR_CALL( updateHashtable(cmr, dense, rowData, processedRows, numProcessedRows, rowHashtable) );
      CMR_CALL( updateHashtable(cmr, dense, columnData, processedColumns, numProcessedColumns, columnHashtable) );
    }
    else
    {
      CMRdbgMsg(8, "No path found!\n");

      /* Find a 1-entry of the rank-1 part. */
      size_t extraRow = SIZE_MAX;
      size_t extraColumn = SIZE_MAX;
      if (CMRelementIsRow(startElement))
      {
        size_t row = CMRelementToRowIndex(startElement);
        for (size_t c = 0; c < numProcessedColumns; ++c)
        {
          size_t column = processedColumns[c];
          if (CMRdensebinmatrixGet(dense, row, column))
          {
            extraRow = row;
            extraColumn = column;
          }
        }
      }
      else
      {
        size_t column = CMRelementToColumnIndex(startElement);
        for (size_t r = 0; r < numProcessedRows; ++r)
        {
          size_t row = processedRows[r];
          if (CMRdensebinmatrixGet(dense, row, column))
          {
            extraRow = row;
            extraColumn = column;
          }
        }
      }
      assert(extraRow < SIZE_MAX);
      assert(extraColumn < SIZE_MAX);

      CMR_SEPA* separation = NULL;
      CMR_CALL( CMRsepaCreate(cmr, numRows, numColumns, &separation) );

      for (size_t row = 0; row < numRows; ++row)
      {
        unsigned char part = (CMRelementIsValid(rowData[row].predecessor) || rowData[row].isSource
          || CMRrowToElement(row) == startElement) ? 1 : 0;
        CMR_ELEMENT original = rowData[row].original;
        if (CMRelementIsRow(original))
          separation->rowsToPart[CMRelementToRowIndex(original)] = part;
        else
          separation->columnsToPart[CMRelementToColumnIndex(original)] = part;
      }
      for (size_t column = 0; column < numColumns; ++column)
      {
        unsigned char part = (CMRelementIsValid(columnData[column].predecessor) || columnData[column].isSource
          || CMRcolumnToElement(column) == startElement) ? 1 : 0;
        CMR_ELEMENT original = columnData[column].original;
        if (CMRelementIsRow(original))
          separation->rowsToPart[CMRelementToRowIndex(original)] = part;
        else
          separation->columnsToPart[CMRelementToColumnIndex(original)] = part;
      }

      size_t extraRow0 = SIZE_MAX;
      size_t extraColumn0 = SIZE_MAX;
      size_t extraRow1 = SIZE_MAX;
      size_t extraColumn1 = SIZE_MAX;

      if (CMRelementIsRow(rowData[extraRow].original))
      {
        assert(CMRelementIsColumn(columnData[extraColumn].original));
        extraRow0 = CMRelementToRowIndex(rowData[extraRow].original);
        extraColumn1 = CMRelementToColumnIndex(columnData[extraColumn].original);
      }
      else
      {
        assert(CMRelementIsColumn(rowData[extraRow].original));
        assert(CMRelementIsRow(columnData[extraColumn].original));
        extraColumn0 = CMRelementToColumnIndex(rowData[extraRow].original);
        extraRow1 = CMRelementToRowIndex(columnData[extraColumn].original);
      }

      CMR_CALL( CMRsepaInitialize(cmr, separation, extraRow0, extraColumn1, extraRow1, extraColumn0, SIZE_MAX, SIZE_MAX,
        SIZE_MAX, SIZE_MAX) );
      CMR_CALL( CMRdecApplySeparation(cmr, dec, separation) );
      CMR_CALL( CMRsepaFree(cmr, &separation) );

      break;
    }
  }

  if (dec->type == CMR_DEC_TWO_SUM)
    CMRdbgMsg(6, "Aborting construction of sequence of nested 3-connected minors due to a 2-separation.\n");
  else
    CMRdbgMsg(6, "Successfully constructed sequence of nested 3-connected minors.\n");

  CMR_CALL( CMRlisthashtableFree(cmr, &columnHashtable) );
  CMR_CALL( CMRlisthashtableFree(cmr, &rowHashtable) );

  CMR_CALL( CMRfreeStackArray(cmr, &processedColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &processedRows) );
  
  CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowData) );
  CMR_CALL( CMRdensebinmatrixFreeStack(cmr, &dense) );
  CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );

  return CMR_OKAY;
}

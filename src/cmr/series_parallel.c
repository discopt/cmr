// #define CMR_DEBUG /* Uncomment to debug this file. */

// TODO: Try to not fill the hashtable in initialScan but just add all elements to the queue.

#include <cmr/series_parallel.h>

#include "env_internal.h"
#include "hashtable.h"
#include "sort.h"
#include "listmatrix.h"

#include <stdint.h>
#include <time.h>

typedef enum
{
  REMOVED = 0,
  ZERO = 1,
  BLOCK = 2,
  OTHER = 3
} ElementType;

typedef struct
{
  long long hashValue;                /**< \brief Hash value of this element. */
  CMR_LISTHASHTABLE_ENTRY hashEntry;  /**< \brief Entry in row or column hashtable. */
  size_t distance;                    /**< \brief Distance in breadth-first search. */
  size_t predecessor;                 /**< \brief Index of predecessor element in breadth-first search. */
  bool inQueue;                       /**< \brief Whether this element is queued. */
  char lastBFS;                       /**< \brief Last breadth-first search that found this node.
                                       **< Is 0 initially, positive for search runs, -1 if marked and -2 for SP-reduced
                                       **< element. */
  bool specialBFS;                    /**< \brief Whether this is a special node in breadth-first search. */
} ElementData;


CMR_ERROR CMRspInitStatistics(CMR_SP_STATISTICS* stats)
{
  assert(stats);

  stats->totalCount = 0;
  stats->totalTime = 0.0;
  stats->reduceCount = 0;
  stats->reduceTime = 0.0;
  stats->nonbinaryCount = 0;
  stats->nonbinaryTime = 0.0;
  stats->wheelCount = 0;
  stats->wheelTime = 0.0;

  return CMR_OKAY;
}

CMR_ERROR CMRspPrintStatistics(FILE* stream, CMR_SP_STATISTICS* stats)
{
  assert(stream);
  assert(stats);

  fprintf(stream, "Series-parallel computations (count / time):\n");
  fprintf(stream, "Search for reductions:          %ld / %f\n", stats->reduceCount, stats->reduceTime);
  fprintf(stream, "Search for wheel matrices:      %ld / %f\n", stats->wheelCount, stats->wheelTime);
  fprintf(stream, "Search for ternary certificate: %ld / %f\n", stats->nonbinaryCount, stats->nonbinaryTime);
  fprintf(stream, "Total:                          %ld / %f\n", stats->totalCount, stats->totalTime);

  return CMR_OKAY;
}

static char seriesParallelStringBuffer[32]; /**< Static buffer for \ref CMRspString. */

char* CMRspReductionString(CMR_SP_REDUCTION reduction, char* buffer)
{
  if (!buffer)
    buffer = seriesParallelStringBuffer;

  if (reduction.element > 0)
  {
    if (reduction.mate > 0)
      sprintf(buffer, "c%d copy of c%d", reduction.element, reduction.mate);
    else if (reduction.mate < 0)
      sprintf(buffer, "c%d unit at r%d", reduction.element, -reduction.mate);
    else
      sprintf(buffer, "c%d zero", reduction.element);
    return buffer;
  }
  else if (reduction.element < 0)
  {
    if (reduction.mate > 0)
      sprintf(buffer, "r%d unit at c%d", -reduction.element, reduction.mate);
    else if (reduction.mate < 0)
      sprintf(buffer, "r%d copy of r%d", -reduction.element, -reduction.mate);
    else
      sprintf(buffer, "r%d zero", -reduction.element);
    return buffer;
  }
  else
    return "<invalid series-parallel reduction>";
}

/**
 * \brief Removes the given nonzero from the linked-list representation.
 */

static inline
void unlinkNonzero(
  ListMatrixNonzero* nonzero /**< Nonzero to be removed from the linked lists. */
)
{
  assert(nonzero);

  CMRdbgMsg(4, "Removing r%d,c%d from linked list.\n", nonzero->row+1, nonzero->column+1);
  nonzero->above->below = nonzero->below;
  nonzero->below->above = nonzero->above;
  nonzero->left->right = nonzero->right;
  nonzero->right->left = nonzero->left;
}

/**
 * \brief Returns the smallest power of 2 at least as large as \p x.
 */

static
size_t nextPower2(size_t x)
{
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return x + 1;
}

/**
 * \brief Allocates and initializes element data (on the stack).
 */

static
CMR_ERROR createElementData(
  CMR* cmr,                   /**< \ref CMR environment. */
  ElementData** pelementData, /**< Pointer for storing the element data. */
  size_t size                 /**< Number of elements. */
)
{
  assert(cmr);
  assert(pelementData);

  CMR_CALL( CMRallocStackArray(cmr, pelementData, size) );
  ElementData* elementData = *pelementData;
  for (size_t e = 0; e < size; ++e)
  {
    elementData[e].hashValue = 0;
    elementData[e].hashEntry = SIZE_MAX;
    elementData[e].inQueue = false;
    elementData[e].lastBFS = 0;
  }

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
    CMRdbgMsg(2, "Entry %d has hash %ld.\n", e, h);
    h = projectSignedHash(3 * h);
  }

  return CMR_OKAY;
}

/**
 * \brief Scan the matrix to compute the number of nonzeros and the hash of each row and each column.
 */

static
CMR_ERROR calcNonzeroCountHashFromMatrix(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,       /**< Matrix. */
  ListMatrix* listmatrix,   /**< List matrix representation. */
  ElementData* rowData,     /**< Other row element data. */
  ElementData* columnData,  /**< Other column element data. */
  long long* hashVector     /**< Hash vector. */
)
{
  assert(cmr);
  assert(matrix);
  assert(rowData);
  assert(columnData);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      char value = matrix->entryValues[e];
      assert(value == 1 || value == -1);

      /* Update row data. */
      listmatrix->rowElements[row].numNonzeros++;
      long long newHash = projectSignedHash(rowData[row].hashValue + value * hashVector[column]);
      rowData[row].hashValue  = newHash;

      /* Update column data. */
      listmatrix->columnElements[column].numNonzeros++;
      newHash = projectSignedHash(columnData[column].hashValue + value * hashVector[row]);
      columnData[column].hashValue = newHash;
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Scan the matrix to compute the number of nonzeros and the hash of each row and each column.
 */

static
CMR_ERROR calcBinaryHashFromListMatrix(
  CMR* cmr,                 /**< \ref CMR environment. */
  ListMatrix* listmatrix,   /**< List matrix. */
  ElementData* rowData,     /**< Other row element data. */
  ElementData* columnData,  /**< Other column element data. */
  long long* hashVector     /**< Hash vector. */
)
{
  assert(cmr);
  assert(listmatrix);
  assert(rowData);
  assert(columnData);

  /* Reset hash values. */
  for (ListMatrixNonzero* rowHead = listmatrix->anchor.below; rowHead->row != SIZE_MAX; rowHead = rowHead->below)
    rowData[rowHead->row].hashValue = 0;
  for (ListMatrixNonzero* columnHead = listmatrix->anchor.right; columnHead->column != SIZE_MAX;
    columnHead = columnHead->right)
  {
    columnData[columnHead->column].hashValue = 0;
  }

  for (ListMatrixNonzero* rowHead = listmatrix->anchor.below; rowHead->row != SIZE_MAX; rowHead = rowHead->below)
  {
    for (ListMatrixNonzero* nz = rowHead->right; nz != rowHead; nz = nz->right)
    {
      /* Update row data. */
      long long newHash = projectSignedHash(rowData[nz->row].hashValue + hashVector[nz->column]);
      rowData[nz->row].hashValue  = newHash;

      /* Update column data. */
      newHash = projectSignedHash(columnData[nz->column].hashValue + hashVector[nz->row]);
      columnData[nz->column].hashValue = newHash;
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Scans the matrix initially in order to add all rows or columns either to the queue or to the hashtable.
 */

static
CMR_ERROR initializeQueueHashtableFromMatrix(
  CMR* cmr,                               /**< \ref CMR environment. */
  CMR_LISTHASHTABLE* hashtable,           /**< Row or column hashtable. */
  ListMatrixElement* listmatrixElements,  /**< Row or column list matrix elements. */
  ElementData* data,                      /**< Other row/column data. */
  size_t sizeData,                        /**< Length of \p ListMatrixElements and \p data. */
  CMR_ELEMENT* queue,                     /**< Queue. */
  size_t* pqueueEnd,                      /**< Pointer to end of queue. */
  bool isRow                              /**< Whether we are deadling with rows. */
)
{
  assert(cmr);
  assert(hashtable || sizeData == 0);

  for (size_t i = 0; i < sizeData; ++i)
  {
    CMRdbgMsg(2, "%s %d has %d nonzeros.\n", isRow ? "Row" : "Column", i, listmatrixElements[i].numNonzeros);

    /* Check if it qualifies for addition to the hashtable. */
    if (listmatrixElements[i].numNonzeros > 1)
    {
      CMR_LISTHASHTABLE_ENTRY entry = CMRlisthashtableFindFirst(hashtable, llabs(data[i].hashValue));
      CMRdbgMsg(2, "Search for hash %ld of %s %d yields entry %d.\n", data[i].hashValue, isRow ? "row" : "column", i, entry);
      if (entry == SIZE_MAX)
      {
        CMR_CALL( CMRlisthashtableInsert(cmr, hashtable, llabs(data[i].hashValue), i, &data[i].hashEntry) );
        continue;
      }
    }

    /* If it was not added to the hashtable, we add it to the queue. */
    queue[*pqueueEnd] = isRow ? CMRrowToElement(i) : CMRcolumnToElement(i);
    data[i].hashEntry = SIZE_MAX;
    data[i].inQueue = true;
    (*pqueueEnd)++;
  }

  return CMR_OKAY;
}

/**
 * \brief Adds all rows or columns of the list representation of the matrix either to the queue or to the hashtable.
 */

static
CMR_ERROR initializeQueueHashtableFromListMatrix(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_LISTHASHTABLE* hashtable, /**< Row or column hashtable. */
  ListMatrix* listmatrix,       /**< List matrix. */
  ElementData* data,            /**< Other row/column data array. */
  CMR_ELEMENT* queue,           /**< Queue. */
  size_t* pqueueEnd,            /**< Pointer to end of queue. */
  bool isRow                    /**< Whether we are deadling with rows. */
)
{
  assert(cmr);

  CMRdbgMsg(2, "Initializing queue and hashtable from list representation. Inspecting %s.\n",
    isRow ? "rows" : "columns");

  ListMatrixNonzero* anchor = &listmatrix->anchor;

  for (ListMatrixNonzero* head = isRow ? anchor->below : anchor->right; head != anchor;
    head = (isRow ? head->below : head->right))
  {
    size_t i = isRow ? head->row : head->column;
    CMRdbgMsg(4, "%s %d has %d nonzeros.\n", isRow ? "Row" : "Column", i,
      isRow ? listmatrix->rowElements[i].numNonzeros : listmatrix->columnElements[i].numNonzeros);

    /* Check if it qualifies for addition to the hashtable. */
    CMR_LISTHASHTABLE_ENTRY entry = CMRlisthashtableFindFirst(hashtable, llabs(data[i].hashValue));
    CMRdbgMsg(6, "Search for hash %ld of %s %d yields entry %d.\n", data[i].hashValue, isRow ? "row" : "column", i, entry);
    if (entry == SIZE_MAX)
      CMR_CALL( CMRlisthashtableInsert(cmr, hashtable, llabs(data[i].hashValue), i, &data[i].hashEntry) );
    else
    {
      /* If it was not added to the hashtable, we add it to the queue. */
      queue[*pqueueEnd] = isRow ? CMRrowToElement(i) : CMRcolumnToElement(i);
      data[i].hashEntry = SIZE_MAX;
      data[i].inQueue = true;
      (*pqueueEnd)++;
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Checks whether the row/column at \p index is a (negated) copy of a row/column stored in the hashtable.
 *
 * Compares the actual vectors by traversing the linked-list representation.
 */

static
size_t findCopy(
  ListMatrixElement* listData,  /**< Row/column list data. */
  ElementData* data,            /**< Other row/column data. */
  CMR_LISTHASHTABLE* hashtable, /**< Row/column hashtable. */
  size_t index,                 /**< Index in \p data. */
  bool isRow,                   /**< Whether we are dealing with rows. */
  bool support                  /**< Whether only the support of the vector matters. */
)
{
  CMR_LISTHASHTABLE_HASH hash = llabs(data[index].hashValue);
  CMRdbgMsg(4, "Processing %s %d with a collision (hash value %ld).\n", isRow ? "row" : "column", index,
    data[index].hashValue);
  for (CMR_LISTHASHTABLE_ENTRY entry = CMRlisthashtableFindFirst(hashtable, hash);
    entry != SIZE_MAX; entry = CMRlisthashtableFindNext(hashtable, hash, entry))
  {
    size_t collisionIndex = CMRlisthashtableValue(hashtable, entry);
    CMRdbgMsg(8, "%s %d has the same hash value. Comparing...\n", isRow ? "Row" : "Column", collisionIndex);
    bool equal = true;
    bool negated = true;
    if (isRow)
    {
      ListMatrixNonzero* nz1 = listData[index].head.right;
      ListMatrixNonzero* nz2 = listData[collisionIndex].head.right;
      while (equal || negated || support)
      {
        if (nz1->column != nz2->column)
        {
          equal = false;
          negated = false;
          support = false;
          break;
        }
        if (nz1->column == SIZE_MAX)
          break;
        if (nz1->value == nz2->value)
          negated = false;
        else
          equal = false;
        nz1 = nz1->right;
        nz2 = nz2->right;
      }
    }
    else
    {
      ListMatrixNonzero* nz1 = listData[index].head.below;
      ListMatrixNonzero* nz2 = listData[collisionIndex].head.below;
      while (equal || negated || support)
      {
        if (nz1->row != nz2->row)
        {
          equal = false;
          negated = false;
          support = false;
          break;
        }
        if (nz1->row == SIZE_MAX)
          break;
        if (nz1->value == nz2->value)
          negated = false;
        else
          equal = false;
        nz1 = nz1->below;
        nz2 = nz2->below;
      }
    }

    if (equal || negated || support)
      return collisionIndex;
  }

  return SIZE_MAX;
}

/**
 * \brief Processes the deletion of a nonzero from the linked-list representation.
 */

static
CMR_ERROR processNonzero(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_LISTHASHTABLE* hashtable,     /**< Row/column hashtable. */
  long long hashChange,             /**< Modification of the hash value. */
  size_t index,                     /**< Index of row/column. */
  ListMatrixElement* indexListData, /**< Row/column list data. */
  ElementData* indexData,           /**< Other row/column data. */
  CMR_ELEMENT* queue,               /**< Queue. */
  size_t* pqueueEnd,                /**< Pointer to end of queue. */
  size_t queueMemory,               /**< Memory allocated for queue. */
  bool isRow                        /**< Whether we are dealing with rows. */
)
{
  assert(cmr);
  assert(indexData);
  assert(hashtable);
  assert(queue);

  indexListData->numNonzeros--;
  long long newHash = projectSignedHash(indexData->hashValue + hashChange);
  CMRdbgMsg(4, "Processing nonzero. Old hash is %ld, change is %ld, new hash is %ld.\n", indexData->hashValue,
    hashChange, newHash);
  indexData->hashValue = newHash;

  /* Add to queue if necessary. */
  if (!indexData->inQueue)
  {
    queue[*pqueueEnd % queueMemory] = isRow ? CMRrowToElement(index) : CMRcolumnToElement(index);
    indexData->inQueue = true;
    (*pqueueEnd)++;
  }

  /* Remove from hashtable if necessary. */
  if (indexData->hashEntry != SIZE_MAX)
  {
    CMR_CALL( CMRlisthashtableRemove(cmr, hashtable, indexData->hashEntry) );
    indexData->hashEntry = SIZE_MAX;
  }

  return CMR_OKAY;
}

/**
 * \brief Carry out the actual reduction algorithm after all data structures are initialized.
 */

static
CMR_ERROR reduceListMatrix(
  CMR* cmr,                           /**< \ref CMR environment. */
  ListMatrix* listmatrix,             /**< List matrix. */
  ElementData* rowData,               /**< Row data. */
  ElementData* columnData,            /**< Column data. */
  CMR_LISTHASHTABLE* rowHashtable,    /**< Row hashtable. */
  CMR_LISTHASHTABLE* columnHashtable, /**< Column hashtable. */
  long long* entryToHash,             /**< Pre-computed hash values of vector entries. */
  CMR_ELEMENT* queue,                 /**< Queue. */
  size_t* pqueueStart,                /**< Pointer to start of queue. */
  size_t* pqueueEnd,                  /**< Pointer to end of queue. */
  size_t queueMemory,                 /**< Memory allocated for queue. */
  CMR_SP_REDUCTION* reductions,       /**< Array for storing the SP-reductions. Must be sufficiently large, e.g., number
                                       **< of rows + number of columns. */
  size_t* pnumReductions,             /**< Pointer for storing the number of SP-reductions. */
  size_t* pnumRowReductions,          /**< Pointer for storing the number of row reductions (may be \c NULL). */
  size_t* pnumColumnReductions        /**< Pointer for storing the number of column reductions (may be \c NULL). */
)
{
  assert(cmr);
  assert(listmatrix);

  while (*pqueueEnd > *pqueueStart)
  {
#if defined(CMR_DEBUG)
    CMRdbgMsg(0, "\n    Status:\n");
    for (ListMatrixNonzero* nz = listmatrix->anchor.below; nz != &listmatrix->anchor; nz = nz->below)
    {
      size_t row = nz->row;
      CMRdbgMsg(6, "Row r%d: %d nonzeros, hashed = %s, hash = %ld", row+1, listmatrix->rowElements[row].numNonzeros,
        rowData[row].hashEntry == SIZE_MAX ? "NO" : "YES", rowData[row].hashValue);
      if (rowData[row].hashEntry != SIZE_MAX)
      {
        CMRdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", rowData[row].hashEntry,
          CMRlisthashtableHash(rowHashtable, rowData[row].hashEntry),
          CMRlisthashtableValue(rowHashtable, rowData[row].hashEntry));
      }
      CMRdbgMsg(0, "\n");
    }
    for (ListMatrixNonzero* nz = listmatrix->anchor.right; nz != &listmatrix->anchor; nz = nz->right)
    {
      size_t column = nz->column;
      CMRdbgMsg(6, "Column c%d: %d nonzeros, hashed = %s, hash = %ld", column+1,
        listmatrix->columnElements[column].numNonzeros, columnData[column].hashEntry == SIZE_MAX ? "NO" : "YES",
        columnData[column].hashValue);
      if (columnData[column].hashEntry != SIZE_MAX)
      {
        CMRdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", columnData[column].hashEntry,
          CMRlisthashtableHash(columnHashtable, columnData[column].hashEntry),
          CMRlisthashtableValue(columnHashtable, columnData[column].hashEntry));
      }
      CMRdbgMsg(0, "\n");
    }

    for (size_t q = *pqueueStart; q < *pqueueEnd; ++q)
    {
      CMR_ELEMENT e = queue[q % queueMemory];
      CMRdbgMsg(6, "Queue #%d @ %d: %s %d\n", q, q % queueMemory,
        CMRelementIsRow(e) ? "row" : "column", CMRelementIsRow(e) ? CMRelementToRowIndex(e) : CMRelementToColumnIndex(e));
    }
    CMRdbgMsg(0, "\n");
#endif /* CMR_DEBUG */

    CMR_ELEMENT element = queue[(*pqueueStart) % queueMemory];
    ++(*pqueueStart);

    CMRdbgMsg(2, "Top element is %s %d with %d nonzeros.\n", CMRelementIsRow(element) ? "row" : "column",
      CMRelementIsRow(element) ? CMRelementToRowIndex(element) : CMRelementToColumnIndex(element),
      CMRelementIsRow(element) ? listmatrix->rowElements[CMRelementToRowIndex(element)].numNonzeros :
      listmatrix->columnElements[CMRelementToColumnIndex(element)].numNonzeros);

    if (CMRelementIsRow(element))
    {
      /* We consider a row. */

      size_t row1 = CMRelementToRowIndex(element);
      if (listmatrix->rowElements[row1].numNonzeros > 1)
      {
        rowData[row1].inQueue = false;
        size_t row2 = findCopy(listmatrix->rowElements, rowData, rowHashtable, row1, true, false);

        if (row2 == SIZE_MAX)
        {
          CMRdbgMsg(6, "No parallel row found. Inserting row %d.\n", row1);
          CMR_CALL( CMRlisthashtableInsert(cmr, rowHashtable, llabs(rowData[row1].hashValue), row1, &rowData[row1].hashEntry) );
        }
        else
        {
          CMRdbgMsg(6, "Row %d is parallel.\n", row2);

          /* We found a parallel row. */
          reductions[*pnumReductions].element = CMRrowToElement(row1);
          reductions[*pnumReductions].mate = CMRrowToElement(row2);
          (*pnumReductions)++;
          (*pnumRowReductions)++;

          for (ListMatrixNonzero* entry = listmatrix->rowElements[row1].head.right; entry->column != SIZE_MAX;
            entry = entry->right)
          {
            CMRdbgMsg(8, "Processing nonzero at column %d.\n", entry->column);

            unlinkNonzero(entry);
            CMR_CALL( processNonzero(cmr, columnHashtable, -entryToHash[entry->row] * entry->value, entry->column,
              &listmatrix->columnElements[entry->column], &columnData[entry->column], queue, pqueueEnd, queueMemory, false) );
          }
          listmatrix->rowElements[row1].numNonzeros = 0;
          rowData[row1].lastBFS = -2;
          listmatrix->rowElements[row1].head.above->below = listmatrix->rowElements[row1].head.below;
          listmatrix->rowElements[row1].head.below->above = listmatrix->rowElements[row1].head.above;

          assert(listmatrix->rowElements[row1].head.left == &listmatrix->rowElements[row1].head);
          assert(listmatrix->rowElements[row1].head.right == &listmatrix->rowElements[row1].head);
        }
      }
      else
      {
        /* Zero or unit row vector. */
        CMRdbgMsg(4, "Processing %s row %d.\n", listmatrix->rowElements[row1].numNonzeros == 0 ? "zero" : "unit", row1);

        rowData[row1].inQueue = false;
        if (listmatrix->rowElements[row1].numNonzeros)
        {
          ListMatrixNonzero* entry = listmatrix->rowElements[row1].head.right;
          size_t column = entry->column;

          CMRdbgMsg(4, "Processing unit row %d with 1 in column %d.\n", row1, column);

          unlinkNonzero(entry);
          listmatrix->rowElements[row1].numNonzeros--;
          CMR_CALL( processNonzero(cmr, columnHashtable, -entryToHash[entry->row] * entry->value, column,
            &listmatrix->columnElements[column], &columnData[column], queue, pqueueEnd, queueMemory, false) );
          reductions[*pnumReductions].mate = CMRcolumnToElement(column);
        }
        else
          reductions[*pnumReductions].mate = 0;
        reductions[*pnumReductions].element = element;
        (*pnumReductions)++;
        (*pnumRowReductions)++;
        rowData[row1].lastBFS = -2;
        listmatrix->rowElements[row1].head.above->below = listmatrix->rowElements[row1].head.below;
        listmatrix->rowElements[row1].head.below->above = listmatrix->rowElements[row1].head.above;
      }
    }
    else
    {
      /* We consider a column. */

      size_t column1 = CMRelementToColumnIndex(element);
      if (listmatrix->columnElements[column1].numNonzeros > 1)
      {
        columnData[column1].inQueue = false;
        size_t column2 = findCopy(listmatrix->columnElements, columnData, columnHashtable, column1, false, false);

        if (column2 == SIZE_MAX)
        {
          CMRdbgMsg(6, "No parallel column found. Inserting column %d.\n", column1);
          CMR_CALL( CMRlisthashtableInsert(cmr, columnHashtable, llabs(columnData[column1].hashValue), column1,
            &columnData[column1].hashEntry) );
        }
        else
        {
          CMRdbgMsg(6, "Column %d is parallel.\n", column2);

          /* We found a parallel column. */
          reductions[*pnumReductions].element = CMRcolumnToElement(column1);
          reductions[*pnumReductions].mate = CMRcolumnToElement(column2);
          (*pnumReductions)++;
          (*pnumColumnReductions)++;
          columnData[column1].lastBFS = -2;

          for (ListMatrixNonzero* entry = listmatrix->columnElements[column1].head.below; entry->row != SIZE_MAX;
            entry = entry->below)
          {
            CMRdbgMsg(8, "Processing nonzero at row %d.\n", entry->row);

            unlinkNonzero(entry);
            CMR_CALL( processNonzero(cmr, rowHashtable, -entryToHash[entry->column] * entry->value, entry->row,
              &listmatrix->rowElements[entry->row], &rowData[entry->row], queue, pqueueEnd, queueMemory, true) );
          }
          listmatrix->columnElements[column1].numNonzeros = 0;
          listmatrix->columnElements[column1].head.left->right = listmatrix->columnElements[column1].head.right;
          listmatrix->columnElements[column1].head.right->left = listmatrix->columnElements[column1].head.left;

          assert(listmatrix->columnElements[column1].head.above == &listmatrix->columnElements[column1].head);
          assert(listmatrix->columnElements[column1].head.below == &listmatrix->columnElements[column1].head);
        }
      }
      else
      {
        /* Zero or unit column vector. */
        CMRdbgMsg(4, "Processing %s column %d.\n",
          listmatrix->columnElements[column1].numNonzeros == 0 ? "zero" : "unit", column1);

        columnData[column1].inQueue = false;
        if (listmatrix->columnElements[column1].numNonzeros)
        {
          ListMatrixNonzero* entry = listmatrix->columnElements[column1].head.below;
          size_t row = entry->row;

          CMRdbgMsg(4, "Processing unit column %d with 1 in row %d.\n", column1, row);

          unlinkNonzero(entry);
          listmatrix->columnElements[column1].numNonzeros--;
          CMR_CALL( processNonzero(cmr, rowHashtable, -entryToHash[entry->column] * entry->value, row,
            &listmatrix->rowElements[row], &rowData[row], queue, pqueueEnd, queueMemory, true) );
          reductions[*pnumReductions].mate = CMRrowToElement(row);
        }
        else
          reductions[*pnumReductions].mate = 0;
        reductions[*pnumReductions].element = element;
        (*pnumReductions)++;
        (*pnumColumnReductions)++;
        columnData[column1].lastBFS = -2;
        listmatrix->columnElements[column1].head.left->right = listmatrix->columnElements[column1].head.right;
        listmatrix->columnElements[column1].head.right->left = listmatrix->columnElements[column1].head.left;
      }
    }
  }

  return CMR_OKAY;
}


/**
 * \brief Try to find an \f$ M_2 \f$-submatrix for a ternary SP-reduced matrix in list representation.
 */

static
CMR_ERROR extractNonbinarySubmatrix(
  CMR* cmr,                           /**< \ref CMR environment. */
  ListMatrix* listmatrix,             /**< List matrix. */
  ElementData* rowData,               /**< Row data. */
  ElementData* columnData,            /**< Column data. */
  CMR_LISTHASHTABLE* rowHashtable,    /**< Row hashtable. */
  CMR_LISTHASHTABLE* columnHashtable, /**< Column hashtable. */
  long long* entryToHash,             /**< Pre-computed hash values of vector entries. */
  CMR_ELEMENT* queue,                 /**< Queue. */
  size_t* pqueueStart,                /**< Pointer to start of queue. */
  size_t* pqueueEnd,                  /**< Pointer to end of queue. */
  size_t queueMemory,                 /**< Memory allocated for queue. */
  CMR_SUBMAT** pviolatorSubmatrix     /**< Pointer for storing the submatrix if an \f$ M_2 \f$-submatrix is found. */
)
{
  assert(cmr);
  assert(listmatrix);
  assert(rowData);
  assert(columnData);
  assert(rowHashtable);
  assert(columnHashtable);
  assert(entryToHash);
  assert(queue);
  assert(pviolatorSubmatrix);

  CMRdbgMsg(2, "Searching for M_2-submatrix.\n");

  *pviolatorSubmatrix = NULL;
  while (*pqueueEnd > *pqueueStart)
  {
#if defined(CMR_DEBUG)
    CMRdbgMsg(0, "\n");
    if (anchor)
    {
      CMRdbgMsg(4, "Status:\n");
      for (ListMatrixNonzero* nz = anchor->below; nz != anchor; nz = nz->below)
      {
        size_t row = nz->row;
        CMRdbgMsg(6, "Row r%d: %d nonzeros, hashed = %s, hash = %ld", row+1, rowListData[row].numNonzeros,
          rowData[row].hashEntry == SIZE_MAX ? "NO" : "YES", rowData[row].hashValue);
        if (rowData[row].hashEntry != SIZE_MAX)
        {
          CMRdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", rowData[row].hashEntry,
            CMRlisthashtableHash(rowHashtable, rowData[row].hashEntry),
            CMRlisthashtableValue(rowHashtable, rowData[row].hashEntry));
        }
        CMRdbgMsg(0, "\n");
      }
      for (ListMatrixNonzero* nz = anchor->right; nz != anchor; nz = nz->right)
      {
        size_t column = nz->column;
        CMRdbgMsg(6, "Column c%d: %d nonzeros, hashed = %s, hash = %ld", column+1, columnListData[column].numNonzeros,
          columnData[column].hashEntry == SIZE_MAX ? "NO" : "YES", columnData[column].hashValue);
        if (columnData[column].hashEntry != SIZE_MAX)
        {
          CMRdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", columnData[column].hashEntry,
            CMRlisthashtableHash(columnHashtable, columnData[column].hashEntry),
            CMRlisthashtableValue(columnHashtable, columnData[column].hashEntry));
        }
        CMRdbgMsg(0, "\n");
      }
    }
    for (size_t q = *pqueueStart; q < *pqueueEnd; ++q)
    {
      CMR_ELEMENT e = queue[q % queueMemory];
      CMRdbgMsg(6, "Queue #%d @ %d: %s\n", q, q % queueMemory, CMRelementString(e, NULL));
    }
    CMRdbgMsg(0, "\n");
#endif /* CMR_DEBUG */

    CMR_ELEMENT element = queue[(*pqueueStart) % queueMemory];
    ++(*pqueueStart);

    CMRdbgMsg(2, "Top element is %s %d with %d nonzeros.\n", CMRelementIsRow(element) ? "row" : "column",
      CMRelementIsRow(element) ? CMRelementToRowIndex(element) : CMRelementToColumnIndex(element),
      CMRelementIsRow(element) ? rowListData[CMRelementToRowIndex(element)].numNonzeros :
      columnListData[CMRelementToColumnIndex(element)].numNonzeros);

    if (CMRelementIsRow(element))
    {
      /* We consider a row. */

      size_t row1 = CMRelementToRowIndex(element);
      assert(listmatrix->rowElements[row1].numNonzeros > 1);
      rowData[row1].inQueue = false;
      size_t row2 = findCopy(listmatrix->rowElements, rowData, rowHashtable, row1, true, true);

      if (row2 == SIZE_MAX)
      {
        CMRdbgMsg(6, "No parallel row found. Inserting row %d.\n", row1);
        CMR_CALL( CMRlisthashtableInsert(cmr, rowHashtable, llabs(rowData[row1].hashValue), row1, &rowData[row1].hashEntry) );
      }
      else
      {
        CMRdbgMsg(6, "Row %d is parallel.\n", row2);

        /* We found a row copy. */
        size_t column1 = SIZE_MAX;
        size_t column2 = SIZE_MAX;
        ListMatrixNonzero* nz1 = listmatrix->rowElements[row1].head.right;
        ListMatrixNonzero* nz2 = listmatrix->rowElements[row2].head.right;
        while (column1 == SIZE_MAX || column2 == SIZE_MAX)
        {
          if (nz1->value == nz2->value)
          {
            if (column1 == SIZE_MAX)
              column1 = nz1->column;
          }
          else
          {
            if (column2 == SIZE_MAX)
              column2 = nz1->column;
          }
          nz1 = nz1->right;
          nz2 = nz2->right;
        }

        CMR_CALL( CMRsubmatCreate(cmr, 2, 2, pviolatorSubmatrix) );
        CMR_SUBMAT* violatorSubmatrix = *pviolatorSubmatrix;
        violatorSubmatrix->rows[0] = row1;
        violatorSubmatrix->rows[1] = row2;
        violatorSubmatrix->columns[0] = column1;
        violatorSubmatrix->columns[1] = column2;
        break;
      }
    }
    else
    {
      /* We consider a column. */

      size_t column1 = CMRelementToColumnIndex(element);
      assert(listmatrix->columnElements[column1].numNonzeros > 1);
      columnData[column1].inQueue = false;
      size_t column2 = findCopy(listmatrix->columnElements, columnData, columnHashtable, column1, false, true);

      if (column2 == SIZE_MAX)
      {
        CMRdbgMsg(6, "No parallel column found. Inserting column %d.\n", column1);
        CMR_CALL( CMRlisthashtableInsert(cmr, columnHashtable, llabs(columnData[column1].hashValue), column1,
          &columnData[column1].hashEntry) );
      }
      else
      {
        CMRdbgMsg(6, "Column %d is parallel.\n", column2);

        /* We found a column copy. */
        size_t row1 = SIZE_MAX;
        size_t row2 = SIZE_MAX;
        ListMatrixNonzero* nz1 = listmatrix->columnElements[column1].head.below;
        ListMatrixNonzero* nz2 = listmatrix->columnElements[column2].head.below;
        while (row1 == SIZE_MAX || row2 == SIZE_MAX)
        {
          if (nz1->value == nz2->value)
          {
            if (row1 == SIZE_MAX)
              row1 = nz1->row;
          }
          else
          {
            if (row2 == SIZE_MAX)
              row2 = nz1->row;
          }
          nz1 = nz1->below;
          nz2 = nz2->below;
        }

        CMR_CALL( CMRsubmatCreate(cmr, 2, 2, pviolatorSubmatrix) );
        CMR_SUBMAT* violatorSubmatrix = *pviolatorSubmatrix;
        violatorSubmatrix->rows[0] = row1;
        violatorSubmatrix->rows[1] = row2;
        violatorSubmatrix->columns[0] = column1;
        violatorSubmatrix->columns[1] = column2;
        break;
      }
    }
  }

  return CMR_OKAY;
}

static
CMR_ERROR breadthFirstSearch(
  CMR* cmr,                           /**< \ref CMR environment. */
  size_t currentBFS,                  /**< Number of this execution of breadth-first search. */
  ListMatrixElement* rowListData,     /**< Row data. */
  ListMatrixElement* columnListData,  /**< Column data. */
  ElementData* rowData,               /**< Row data. */
  ElementData* columnData,            /**< Column data. */
  CMR_ELEMENT* queue,                 /**< Queue. */
  size_t queueMemory,                 /**< Memory for queue. */
  CMR_ELEMENT* sources,               /**< Array of source nodes. */
  size_t numSources,                  /**< Number of source nodes. */
  CMR_ELEMENT* targets,               /**< Array of target nodes. */
  size_t numTargets,                  /**< Number of target nodes. */
  size_t* pfoundTarget,               /**< Pointer for storing the index of the target node found. */
  size_t* pnumEdges                   /**< Pointer for storing the number of traversed edges. */
)
{
  assert(cmr);
  assert(rowData);
  assert(columnData);
  assert(queue);
  assert(queueMemory > 0);

  CMRdbgMsg(6, "BFS #%d for %d sources and %d targets.\n", currentBFS, numSources, numTargets);
  size_t queueStart = 0;
  size_t queueEnd = numSources;
  for (size_t s = 0; s < numSources; ++s)
  {
    CMR_ELEMENT element = sources[s];
    queue[s] = element;
    if (CMRelementIsRow(element))
    {
      size_t row = CMRelementToRowIndex(element);
      rowData[row].distance = 0;
      rowData[row].lastBFS = currentBFS;
    }
    else
    {
      size_t column = CMRelementToColumnIndex(element);
      columnData[column].distance = 0;
      columnData[column].lastBFS = currentBFS;
    }
  }
  for (size_t t = 0; t < numTargets; ++t)
  {
    CMR_ELEMENT element = targets[t];
    if (CMRelementIsRow(element))
      rowData[CMRelementToRowIndex(element)].lastBFS = currentBFS+1;
    else
      columnData[CMRelementToColumnIndex(element)].lastBFS = currentBFS+1;
  }
  if (pnumEdges)
    *pnumEdges = 0;
  while (queueEnd > queueStart)
  {
    CMR_ELEMENT element = queue[queueStart % queueMemory];
    queueStart++;
    
    CMRdbgMsg(8, "Queue: %s", CMRelementString(element, NULL));
    for (size_t q = queueStart; q < queueEnd; ++q)
      CMRdbgMsg(0, ",%s", CMRelementString(queue[q % queueMemory], NULL));
    CMRdbgMsg(0, "\n");

    if (CMRelementIsRow(element))
    {
      size_t row = CMRelementToRowIndex(element);
      for (ListMatrixNonzero* nz = rowListData[row].head.right; nz->column != SIZE_MAX; nz = nz->right)
      {
        /* Skip edge if disabled. */
        if (nz->special)
        {
          CMRdbgMsg(10, "Edge r%d,c%d is disabled.\n", nz->row+1, nz->column+1);
          continue;
        }

        if (pnumEdges)
          (*pnumEdges)++;
        size_t column = nz->column;
        if (columnData[column].lastBFS != currentBFS)
        {
          /* We found a new column node. */
          CMRdbgMsg(10, "Node c%d receives distance %d (last BFS #%d)\n", column+1, rowData[row].distance + 1,
            columnData[column].lastBFS);
          columnData[column].distance = rowData[row].distance + 1;
          columnData[column].predecessor = row;
          queue[queueEnd % queueMemory] = CMRcolumnToElement(column);
          queueEnd++;
          if (columnData[column].lastBFS == currentBFS+1)
          {
            queueStart = queueEnd;
            if (pfoundTarget)
            {
              for (size_t t = 0; t < numTargets; ++t)
              {
                if (targets[t] == CMRcolumnToElement(column))
                {
                  *pfoundTarget = t;
                  break;
                }
              }
            }
            break;
          }
          else
            columnData[column].lastBFS = currentBFS;
        }
        else
        {
          CMRdbgMsg(10, "Node c%d already known.\n", column+1);
        }
      }
    }
    else
    {
      size_t column = CMRelementToColumnIndex(element);
      for (ListMatrixNonzero* nz = columnListData[column].head.below; nz->row != SIZE_MAX; nz = nz->below)
      {
        /* Skip edge if disabled. */
        if (nz->special)
        {
          CMRdbgMsg(10, "Edge r%d,c%d is disabled.\n", nz->row+1, nz->column+1);
          continue;
        }

        if (pnumEdges)
          (*pnumEdges)++;
        size_t row = nz->row;
        if (rowData[row].lastBFS != currentBFS)
        {
          /* We found a new row node. */
          CMRdbgMsg(10, "Node r%d receives distance %d (last BFS #%d)\n", row+1, columnData[column].distance + 1,
            rowData[row].lastBFS);
          rowData[row].distance = columnData[column].distance + 1;
          rowData[row].predecessor = column;
          queue[queueEnd % queueMemory] = CMRrowToElement(row);
          queueEnd++;
          if (rowData[row].lastBFS == currentBFS+1)
          {
            queueStart = queueEnd;
            if (pfoundTarget)
            {
              for (size_t t = 0; t < numTargets; ++t)
              {
                if (targets[t] == CMRrowToElement(row))
                {
                  *pfoundTarget = t;
                  break;
                }
              }
            }
            break;
          }
          else
            rowData[row].lastBFS = currentBFS;
        }
        else
        {
          CMRdbgMsg(10, "Node r%d already known.\n", row+1);
        }
      }
    }
  }

  /* Reset lastBFS for targets. */
  for (size_t t = 0; t < numTargets; ++t)
  {
    CMR_ELEMENT element = targets[t];
    if (CMRelementIsRow(element))
      rowData[CMRelementToRowIndex(element)].lastBFS = 0;
    else
      columnData[CMRelementToColumnIndex(element)].lastBFS = 0;
  }
  if (pnumEdges)
    (*pnumEdges) /= 2;

  return CMR_OKAY;
}

/**
 * \brief Extract remaining submatrix.
 */

static
CMR_ERROR extractRemainingSubmatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Matrix. */
  size_t numRowReductions,        /**< Number of row SP reductions. */
  size_t numColumnReductions,     /**< Number of column SP reductions. */
  ListMatrix* listmatrix,         /**< List matrix. */
  CMR_SUBMAT** preducedSubmatrix  /**< Pointer for storing the reduced submatrix. */
)
{
  assert(cmr);
  assert(matrix);
  assert(listmatrix);
  assert(preducedSubmatrix);

  CMR_CALL( CMRsubmatCreate(cmr, matrix->numRows - numRowReductions, matrix->numColumns - numColumnReductions,
    preducedSubmatrix) );
  CMR_SUBMAT* remainingSubmatrix = *preducedSubmatrix;
  size_t rowSubmatrix = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (listmatrix->rowElements[row].numNonzeros > 0)
      remainingSubmatrix->rows[rowSubmatrix++] = row;
  }
  assert(rowSubmatrix + numRowReductions == matrix->numRows);

  size_t columnSubmatrix = 0;
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    if (listmatrix->columnElements[column].numNonzeros > 0)
      remainingSubmatrix->columns[columnSubmatrix++] = column;
  }
  assert(columnSubmatrix + numColumnReductions == matrix->numColumns);

  return CMR_OKAY;
}

/**
 * \brief Create full matrix as remaining submatrix.
 */

static
CMR_ERROR createFullRemainingMatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Matrix. */
  CMR_SUBMAT** preducedSubmatrix  /**< Pointer for storing the reduced submatrix. */
)
{
  CMR_CALL( CMRsubmatCreate(cmr, matrix->numRows, matrix->numColumns, preducedSubmatrix) );
  CMR_SUBMAT* remainingSubmatrix = *preducedSubmatrix;
  for (size_t row = 0; row < matrix->numRows; ++row)
    remainingSubmatrix->rows[row] = row;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    remainingSubmatrix->columns[column] = column;

  return CMR_OKAY;
}

/**
 * \brief Searches for a wheel-submatrix.
 */

static
CMR_ERROR extractWheelSubmatrix(
  CMR* cmr,                     /**< \ref CMR environment. */
  ListMatrix* listmatrix,       /**< List matrix. */
  ElementData* rowData,         /**< Row element data. */
  ElementData* columnData,      /**< Column element data. */
  CMR_ELEMENT* queue,           /**< Queue memory. */
  size_t queueMemory,           /**< Size of \p queue. */
  size_t numReducedRows,        /**< Number of rows in reduced matrix. */
  size_t numReducedColumns,     /**< Number of columns in reduced matrix. */
  CMR_SUBMAT** pwheelSubmatrix, /**< Pointer for storing the wheel submatrix (may be \c NULL). */
  CMR_SEPA** pseparation        /**< Pointer for storing a 2-separation (may be \c NULL). */
)
{
  assert(cmr);
  assert(listmatrix);

  CMRdbgMsg(2, "Searching for wheel graph representation submatrix.\n");

  size_t numRows = listmatrix->numRows;
  size_t numColumns = listmatrix->numColumns;
  
  size_t currentBFS = 0;
  for (size_t row = 0; row < numRows; ++row)
    rowData[row].specialBFS = false;
  for (size_t column = 0; column < numColumns; ++column)
    columnData[column].specialBFS = false;
  ListMatrixNonzero** nzBlock = NULL; /* Pointers for simultaneously traversing columns of block. */
  CMR_CALL( CMRallocStackArray(cmr, &nzBlock, numColumns) );
  CMR_ELEMENT* sources = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &sources, numRows) );
  CMR_ELEMENT* targets = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &targets, numColumns) );
  size_t numEdges = 0;
  for (size_t row = 0; row < numRows; ++row)
    numEdges += listmatrix->rowElements[row].numNonzeros;
  while (true)
  {
    size_t sourceRow = listmatrix->anchor.below->row;
    rowData[sourceRow].lastBFS = currentBFS;
    rowData[sourceRow].predecessor = SIZE_MAX;
    rowData[sourceRow].distance = 0;
    size_t targetColumn = listmatrix->rowElements[sourceRow].head.right->column;
    listmatrix->rowElements[sourceRow].head.right->special = 1;

    CMRdbgMsg(4, "Searching for a chordless cycle from r%d to c%d.\n", sourceRow+1, targetColumn+1);

    sources[0] = CMRrowToElement(sourceRow);
    targets[0] = CMRcolumnToElement(targetColumn);
    size_t foundTarget = SIZE_MAX;
    currentBFS++;
    CMR_CALL( breadthFirstSearch(cmr, currentBFS, listmatrix->rowElements, listmatrix->columnElements, rowData,
      columnData, queue, queueMemory, sources, 1, targets, 1, &foundTarget, 0) );
    listmatrix->rowElements[sourceRow].head.right->special = 0;
    assert(foundTarget == 0);
    size_t length = columnData[targetColumn].distance + 1;

    CMRdbgMsg(4, "Length of cycle is %d.\n", length);

    if (length > 4)
    {
      /* We found a long chordless cycle. Traverse backwards along path and collect rows/columns. */
      if (pwheelSubmatrix)
      {
        CMR_CALL( CMRsubmatCreate(cmr, length/2, length/2, pwheelSubmatrix) );
        CMR_SUBMAT* wheelSubmatrix = *pwheelSubmatrix;
        size_t column = targetColumn;
        size_t row;
        for (size_t i = 0; i < length/2; ++i)
        {
          row = columnData[column].predecessor;
          wheelSubmatrix->columns[i] = column;
          wheelSubmatrix->rows[i] = row;
          column = rowData[row].predecessor;
        }
      }
      break;
    }

    /* We have found a 2-by-2 matrix with only 1's. */
    size_t row1 = sourceRow;
    size_t row2 = columnData[targetColumn].predecessor;
    CMRdbgMsg(4, "Growing the 2x2 submatrix with 1's is at r%d, r%d, c%d, c%d.\n", row1+1, row2+1,
      targetColumn+1, rowData[row2].predecessor+1);

    /* Go trough the two nonzeros of the two rows simultaneously. */
    ListMatrixNonzero* nz1 = listmatrix->rowElements[row1].head.right;
    ListMatrixNonzero* nz2 = listmatrix->rowElements[row2].head.right;
    size_t numTargets = 0;
    while (nz1->column != SIZE_MAX)
    {
      CMRdbgMsg(6, "nonzeros at column indices %d and %d.\n", nz1->column, nz2->column);
      if (nz1->column < nz2->column)
      {
        nz1 = nz1->right;
      }
      else if (nz1->column > nz2->column)
        nz2 = nz2->right;
      else
      {
        columnData[nz1->column].specialBFS = true;
        nzBlock[numTargets] = listmatrix->columnElements[nz1->column].head.below;
        targets[numTargets++] = CMRcolumnToElement(nz1->column);
        nz1 = nz1->right;
        nz2 = nz2->right;
      }
    }
    CMRdbgMsg(4, "Identified %d target columns.\n", numTargets);
    assert(numTargets >= 2);

    /* Go through the nonzeros of all marked columns simultaneously. */
    size_t maxIndex = 0;
    size_t maxRow = nzBlock[0]->row;
    size_t currentIndex = 1;
    size_t numSources = 0;
    while (maxRow != SIZE_MAX)
    {
      if (currentIndex == maxIndex)
      {
        /* All nonzeros now have the same row. */
        for (size_t j = 0; j < numTargets; ++j)
          nzBlock[j]->special = 1;
        rowData[maxRow].specialBFS = true;
        sources[numSources++] = CMRrowToElement(maxRow);
        ++maxRow;
      }
      while (nzBlock[currentIndex]->row < maxRow)
        nzBlock[currentIndex] = nzBlock[currentIndex]->below;
      if (nzBlock[currentIndex]->row > maxRow)
      {
        maxIndex = currentIndex;
        maxRow = nzBlock[currentIndex]->row;
      }
      currentIndex = (currentIndex + 1) % numTargets;
    }

    CMRdbgMsg(4, "Identified %d source rows.\n", numSources);
    assert(numSources >= 2);

    currentBFS++;
    foundTarget = SIZE_MAX;
    size_t numTraversedEdges = 0;
    CMR_CALL( breadthFirstSearch(cmr, currentBFS, listmatrix->rowElements, listmatrix->columnElements, rowData,
      columnData, queue, queueMemory, sources, numSources, targets, numTargets, &foundTarget, &numTraversedEdges) );

    if (foundTarget < SIZE_MAX && pwheelSubmatrix)
    {
      size_t length = columnData[CMRelementToColumnIndex(targets[foundTarget])].distance + 1;
      targetColumn = CMRelementToColumnIndex(targets[foundTarget]);
      CMRdbgMsg(4, "Length of cycle is %d.\n", length);

      size_t wheelSize = length == 4 ? 3 : length/2;
      CMR_CALL( CMRsubmatCreate(cmr, wheelSize, wheelSize, pwheelSubmatrix) );
      CMR_SUBMAT* wheelSubmatrix = *pwheelSubmatrix;

      if (length > 4)
      {
        /* The cycle is long and we simply traverse it and collect all rows and columns. */
        size_t column = targetColumn;
        size_t row;
        for (size_t i = 0; i < length/2; ++i)
        {
          row = columnData[column].predecessor;
          wheelSubmatrix->columns[i] = column;
          wheelSubmatrix->rows[i] = row;
          column = rowData[row].predecessor;
        }
      }
      else
      {
        /* The cycle is short and we need to add one row and one column from the block of 1's. */

        size_t column1 = targetColumn; /* Belongs to block. */
        size_t row2 = columnData[column1].predecessor; /* Does not belong to block. */
        size_t column2 = rowData[row2].predecessor; /* Does not belong to block. */
        size_t row1 = columnData[column2].predecessor; /* Belongs to block. */
        CMRdbgMsg(4, "Short cycle is induced by r%d,r%d,c%d,c%d.\n", row1+1, row2+1, column1+1, column2+1);

        /* Go through nonzeros of row2 and mark them. */
        for (ListMatrixNonzero* nz = listmatrix->rowElements[row2].head.right; nz->column != SIZE_MAX; nz = nz->right)
          columnData[nz->column].lastBFS = -1;

        /* Find a non-marked source column. */
        size_t column3 = SIZE_MAX;
        for (size_t t = 0; t < numTargets; ++t)
        {
          column3 = CMRelementToColumnIndex(targets[t]);
          if (columnData[column3].lastBFS >= 0)
            break;
        }
        assert(column3 != SIZE_MAX);
        CMRdbgMsg(4, "Adding c%d\n", column3+1);

        /* Go through nonzeros of column 2 and mark them. */
        for (ListMatrixNonzero* nz = listmatrix->columnElements[column2].head.below; nz->row != SIZE_MAX; nz = nz->below)
          rowData[nz->row].lastBFS = -1;

        /* Find a non-marked source column. */
        size_t row3 = SIZE_MAX;
        for (size_t s = 0; s < numSources; ++s)
        {
          row3 = CMRelementToRowIndex(sources[s]);
          if (rowData[row3].lastBFS >= 0)
            break;
        }
        assert(row3 != SIZE_MAX);
        CMRdbgMsg(4, "Adding r%d\n", row3+1);
        wheelSubmatrix->rows[0] = row3;
        wheelSubmatrix->rows[1] = row1;
        wheelSubmatrix->rows[2] = row2;
        wheelSubmatrix->columns[0] = column3;
        wheelSubmatrix->columns[1] = column1;
        wheelSubmatrix->columns[2] = column2;
      }
      break;
    }
    else if (foundTarget < SIZE_MAX)
    {
      /* Wheel found, but the user does not care. */
      break;
    }

    if (pseparation)
    {
      CMRdbgMsg(4, "No path found. Extracting 2-separation.\n");

      CMR_CALL( CMRsepaCreate(cmr, numReducedRows, numReducedColumns, pseparation) );
      CMR_SEPA* sepa = *pseparation;

      /* Collect all reduced reachable rows/columns. */
      size_t reducedRow = 0;
      size_t reducedSourceRow = 0;
      for (size_t row = 0; row < numRows; ++row)
      {
        if (rowData[row].lastBFS == -2)
          continue;

        sepa->rowsToPart[reducedRow] = (rowData[row].lastBFS == currentBFS) ? 0 : 1;
        CMRdbgMsg(6, "Assigning row r%ld = reduced row r%ld to part %d.\n", row+1, reducedRow+1,
          sepa->rowsToPart[reducedRow]);
        if (row == sourceRow)
          reducedSourceRow = reducedRow;
        ++reducedRow;
      }

      /* Collect all reduced columns that are not reachable. */
      size_t reducedColumn = 0;
      size_t reducedTargetColumn = 0;
      for (size_t column = 0; column < numColumns; ++column)
      {
        if (columnData[column].lastBFS == -2)
          continue;

        sepa->columnsToPart[reducedColumn] = (columnData[column].lastBFS == currentBFS) ? 0 : 1;
        CMRdbgMsg(6, "Assigning column c%ld = reduced column c%ld to part %d.\n", column+1, reducedColumn+1,
          sepa->columnsToPart[column]);
        if (column == targetColumn)
          reducedTargetColumn = reducedColumn;
        ++reducedColumn;
      }

      CMRdbgMsg(6, "Extra row r%ld = reduced row r%ld for part 1 and extra column c%ld = reduced column c%ld for part 0.\n",
        sourceRow + 1, reducedSourceRow + 1, targetColumn + 1, reducedTargetColumn + 1);

      CMR_CALL( CMRsepaInitialize(cmr, sepa, SIZE_MAX, SIZE_MAX, reducedSourceRow, reducedTargetColumn, SIZE_MAX,
        SIZE_MAX, SIZE_MAX, SIZE_MAX) );

      break;
    }

    CMRdbgMsg(4, "No path found. Determining smaller part of 2-separation.\n");
    CMRdbgMsg(4, "%d edges were traversed, %d are in rank-1 submatrix, and %d of %d remain.\n", numTraversedEdges,
      numSources * numTargets, numEdges - numTraversedEdges - numSources * numTargets, numEdges);
    size_t numEdgesPartWithColumn = numTraversedEdges + numSources;
    size_t numEdgesPartWithRow = numEdges - numTraversedEdges - (numSources - 1) * numTargets;
    CMRdbgMsg(4, "Part with a column from block has %d nonzeros, part with a row from the block has %d nonzeros.\n",
      numEdgesPartWithColumn, numEdgesPartWithRow);

    /* We restrict our search on the smaller part. */
    if (numEdgesPartWithColumn <= numEdgesPartWithRow)
    {
      numEdges = numEdgesPartWithColumn;

      /* Remove all but one column from the block. */
      for (size_t t = 0; t < numTargets; ++t)
      {
        size_t column = CMRelementToColumnIndex(targets[t]);
        for (ListMatrixNonzero* nz = listmatrix->columnElements[column].head.below; nz->row != SIZE_MAX; nz = nz->below)
        {
          if (t > 0 || !rowData[nz->row].specialBFS)
          {
            listmatrix->rowElements[nz->row].numNonzeros--;
            unlinkNonzero(nz);
          }
          else
            nz->special = 0;
        }
        if (t == 0)
          listmatrix->columnElements[column].numNonzeros = numSources;
        else
        {
          listmatrix->columnElements[column].numNonzeros = 0;
          listmatrix->columnElements[column].head.left->right = listmatrix->columnElements[column].head.right;
          listmatrix->columnElements[column].head.right->left = listmatrix->columnElements[column].head.left;
        }
      }
      
      /* Mark all source rows and the single remaining target column as non-special again. */
      for (size_t s = 0; s < numSources; ++s)
        rowData[CMRelementToRowIndex(sources[s])].specialBFS = false;
      columnData[CMRelementToColumnIndex(targets[0])].specialBFS = false;
    }
    else
    {
      numEdges = numEdgesPartWithRow;

      /* Remove all but one row from the block. */
      for (size_t s = 0; s < numSources; ++s)
      {
        size_t row = CMRelementToRowIndex(sources[s]);
        for (ListMatrixNonzero* nz = listmatrix->rowElements[row].head.right; nz->column != SIZE_MAX; nz = nz->right)
        {
          if (s > 0 || !columnData[nz->column].specialBFS)
          {
            listmatrix->columnElements[nz->column].numNonzeros--;
            unlinkNonzero(nz);
          }
          else
            nz->special = 0;
        }
        if (s == 0)
          listmatrix->rowElements[row].numNonzeros = numTargets;
        else
        {
          listmatrix->rowElements[row].numNonzeros = 0;
          listmatrix->rowElements[row].head.left->right = listmatrix->rowElements[row].head.right;
          listmatrix->rowElements[row].head.right->left = listmatrix->rowElements[row].head.left;
        }
      }
      
      /* Mark all target columns and the single remaining source row as non-special again. */
      for (size_t t = 0; t < numTargets; ++t)
        columnData[CMRelementToColumnIndex(targets[t])].specialBFS = false;
      rowData[CMRelementToRowIndex(sources[0])].specialBFS = false;
    }
    
    CMRdbgMsg(0, "!!! Recursing.\n");
  }
  CMR_CALL( CMRfreeStackArray(cmr, &targets) );
  CMR_CALL( CMRfreeStackArray(cmr, &sources) );
  CMR_CALL( CMRfreeStackArray(cmr, &nzBlock) );

  return CMR_OKAY;
}

static
CMR_ERROR decomposeBinarySeriesParallel(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Sparse char matrix. */
  CMR_SP_REDUCTION* reductions,     /**< Array for storing the SP-reductions. Must have capacity at least number of
                                     **< rows + number of columns. */
  size_t* pnumReductions,           /**< Pointer for storing the number of SP-reductions. */
  CMR_SUBMAT** preducedSubmatrix,   /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix,  /**< Pointer for storing a wheel-submatrix (may be \c NULL). */
  CMR_SEPA** pseparation,           /**< Pointer for storing a 2-separation (may be \c NULL). */
  CMR_SP_STATISTICS* stats          /**< Pointer to statistics (may be \c NULL). */
)
{
  assert(cmr);
  assert(matrix);
  assert(reductions && pnumReductions);

  CMRdbgMsg(0, "decomposeBinarySeriesParallel for a %dx%d matrix with %d nonzeros.\n", matrix->numRows,
    matrix->numColumns, matrix->numNonzeros);

  clock_t totalClock = 0;
  clock_t reduceClock = 0;
  if (stats)
  {
    reduceClock = clock();
    totalClock = clock();
  }

  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;
  
  /* Create list matrix to use numNonzeros. */
  ListMatrix* listmatrix = NULL;
  CMR_CALL( CMRlistmatrixAlloc(cmr, numRows, numColumns, matrix->numNonzeros, &listmatrix) );
  for (size_t row = 0; row < numRows; ++row)
    listmatrix->rowElements[row].numNonzeros = 0;
  for (size_t column = 0; column < numColumns; ++column)
    listmatrix->columnElements[column].numNonzeros = 0;

  /* Initialize element data and hash vector. */
  ElementData* rowData = NULL;
  CMR_CALL( createElementData(cmr, &rowData, matrix->numRows) );
  ElementData* columnData = NULL;
  CMR_CALL( createElementData(cmr, &columnData, matrix->numColumns) );
  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector,
    matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns) );

  /* Scan the matrix to initialize the element data. */
  CMR_CALL( calcNonzeroCountHashFromMatrix(cmr, matrix, listmatrix, rowData, columnData, hashVector) );

  /* Initialize the queue. */
  CMR_ELEMENT* queue = NULL;
  size_t queueStart = 0;
  size_t queueEnd = 0;
  size_t queueMemory = numRows + numColumns;
  CMR_CALL( CMRallocStackArray(cmr, &queue, queueMemory) );

  /* Initialize the hashtables. */
  CMR_LISTHASHTABLE* rowHashtable = NULL;
  if (numRows > 0)
  {
    CMR_CALL( CMRlisthashtableCreate(cmr, &rowHashtable, nextPower2(numRows), numRows) );
    CMR_CALL( initializeQueueHashtableFromMatrix(cmr, rowHashtable, listmatrix->rowElements, rowData, numRows, queue, &queueEnd,
      true) );
  }
  CMR_LISTHASHTABLE* columnHashtable = NULL;
  if (numColumns > 0)
  {
    CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns) );
    CMR_CALL( initializeQueueHashtableFromMatrix(cmr, columnHashtable, listmatrix->columnElements, columnData, numColumns, queue,
      &queueEnd, false) );
  }

  *pnumReductions = 0;
  if (queueEnd > queueStart || (pviolatorSubmatrix && (numRows + numColumns > 0)))
  {
    /* Initialize list matrix representation. */
    CMR_CALL( CMRlistmatrixInitializeFromMatrix(cmr, listmatrix, matrix) );

    /* We now start main loop. */
    size_t numRowReductions = 0;
    size_t numColumnReductions = 0;
    CMR_CALL( reduceListMatrix(cmr, listmatrix, rowData, columnData, rowHashtable, columnHashtable, hashVector,
      queue, &queueStart, &queueEnd, queueMemory, reductions, pnumReductions, &numRowReductions,
      &numColumnReductions) );

    if (stats)
    {
      stats->reduceCount++;
      stats->reduceTime += (clock() - reduceClock) * 1.0 / CLOCKS_PER_SEC;
    }

    /* Extract remaining submatrix. */
    if (preducedSubmatrix)
    {
      CMR_CALL( extractRemainingSubmatrix(cmr, matrix, numRowReductions, numColumnReductions, listmatrix,
        preducedSubmatrix) );
    }

    if ((pviolatorSubmatrix || pseparation) && (*pnumReductions != (matrix->numRows + matrix->numColumns)))
    {
      clock_t wheelClock = 0;
      if (stats)
        wheelClock = clock();

      CMR_CALL( extractWheelSubmatrix(cmr, listmatrix, rowData, columnData, queue, queueMemory,
        numRows - numRowReductions, numColumns - numColumnReductions, pviolatorSubmatrix, pseparation) );

      if (stats)
      {
        stats->wheelCount++;
        stats->wheelTime += (clock() - wheelClock) * 1.0 / CLOCKS_PER_SEC;
      }
    }
  }
  else
  {
    CMRdbgMsg(2, "No series/parallel element found.\n");

    if (preducedSubmatrix)
      CMR_CALL( createFullRemainingMatrix(cmr, matrix, preducedSubmatrix) );
  }

  CMR_CALL( CMRlisthashtableFree(cmr, &columnHashtable) );
  CMR_CALL( CMRlisthashtableFree(cmr, &rowHashtable) );

  CMR_CALL( CMRfreeStackArray(cmr, &queue) );
  CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowData) );
  
  CMR_CALL( CMRlistmatrixFree(cmr, &listmatrix) );

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestBinarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool* pisSeriesParallel,
  CMR_SP_REDUCTION* reductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix,
  CMR_SP_STATISTICS* stats)
{
  assert(cmr);
  assert(matrix);
  assert(reductions || !pnumReductions);
  assert(!reductions || pnumReductions);

  CMR_SP_REDUCTION* localReductions = NULL;
  size_t localNumReductions = 0;
  if (!reductions)
    CMR_CALL( CMRallocStackArray(cmr, &localReductions, matrix->numRows + matrix->numColumns) );

  CMR_CALL( decomposeBinarySeriesParallel(cmr, matrix, reductions ? reductions : localReductions,
    &localNumReductions, preducedSubmatrix, pviolatorSubmatrix, NULL, stats) );

  if (pisSeriesParallel)
    *pisSeriesParallel = (*pnumReductions == matrix->numRows + matrix->numColumns);
  if (reductions)
    *pnumReductions = localNumReductions;
  else
    CMR_CALL( CMRfreeStackArray(cmr, &localReductions) );

  return CMR_OKAY;
}


CMR_ERROR CMRdecomposeBinarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool* pisSeriesParallel,
  CMR_SP_REDUCTION* reductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix,
  CMR_SEPA** pseparation, CMR_SP_STATISTICS* stats)
{
  assert(cmr);
  assert(matrix);
  assert(reductions || !pnumReductions);
  assert(!reductions || pnumReductions);

  CMR_SP_REDUCTION* localReductions = NULL;
  size_t localNumReductions = 0;
  if (!reductions)
    CMR_CALL( CMRallocStackArray(cmr, &localReductions, matrix->numRows + matrix->numColumns) );

  CMR_CALL( decomposeBinarySeriesParallel(cmr, matrix, reductions ? reductions : localReductions,
    &localNumReductions, preducedSubmatrix, pviolatorSubmatrix, pseparation, stats) );

  if (pisSeriesParallel)
    *pisSeriesParallel = (localNumReductions == matrix->numRows + matrix->numColumns);
  if (reductions)
    *pnumReductions = localNumReductions;
  else
    CMR_CALL( CMRfreeStackArray(cmr, &localReductions) );

  return CMR_OKAY;
}

static
CMR_ERROR decomposeTernarySeriesParallel(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Sparse char matrix. */
  CMR_SP_REDUCTION* reductions,     /**< Array for storing the SP-reductions. Must have capacity at least number of
                                     **< rows + number of columns. */
  size_t* pnumReductions,           /**< Pointer for storing the number of SP-reductions. */
  CMR_SUBMAT** preducedSubmatrix,   /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix,  /**< Pointer for storing a wheel-submatrix (may be \c NULL). */
  CMR_SEPA** pseparation,           /**< Pointer for storing a 2-separation (may be \c NULL). */
  CMR_SP_STATISTICS* stats          /**< Pointer to statistics (may be \c NULL). */
)
{
  assert(cmr);
  assert(matrix);
  assert(reductions && pnumReductions);

  CMRdbgMsg(0, "decomposeTernarySeriesParallel for a %dx%d matrix with %d nonzeros.\n", matrix->numRows,
    matrix->numColumns, matrix->numNonzeros);

  clock_t totalClock = 0;
  clock_t reduceClock = 0;
  if (stats)
  {
    reduceClock = clock();
    totalClock = clock();
  }

  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;

  /* Create list matrix to use numNonzeros. */
  ListMatrix* listmatrix = NULL;
  CMR_CALL( CMRlistmatrixAlloc(cmr, numRows, numColumns, matrix->numNonzeros, &listmatrix) );
  for (size_t row = 0; row < numRows; ++row)
    listmatrix->rowElements[row].numNonzeros = 0;
  for (size_t column = 0; column < numColumns; ++column)
    listmatrix->columnElements[column].numNonzeros = 0;

  /* Initialize element data and hash vector. */
  ElementData* rowData = NULL;
  CMR_CALL( createElementData(cmr, &rowData, matrix->numRows) );
  ElementData* columnData = NULL;
  CMR_CALL( createElementData(cmr, &columnData, matrix->numColumns) );
  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector,
    matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns) );

  /* Scan the matrix to initialize the element data. */
  CMR_CALL( calcNonzeroCountHashFromMatrix(cmr, matrix, listmatrix, rowData, columnData, hashVector) );

  /* Initialize the queue. */
  CMR_ELEMENT* queue = NULL;
  size_t queueStart = 0;
  size_t queueEnd = 0;
  size_t queueMemory = numRows + numColumns;
  CMR_CALL( CMRallocStackArray(cmr, &queue, queueMemory) );

  /* Initialize the hashtables. */
  CMR_LISTHASHTABLE* rowHashtable = NULL;
  if (numRows > 0)
  {
    CMR_CALL( CMRlisthashtableCreate(cmr, &rowHashtable, nextPower2(numRows), numRows) );
    CMR_CALL( initializeQueueHashtableFromMatrix(cmr, rowHashtable, listmatrix->rowElements, rowData, numRows, queue, &queueEnd, true) );
  }
  CMR_LISTHASHTABLE* columnHashtable = NULL;
  if (numColumns > 0)
  {
    CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns) );
    CMR_CALL( initializeQueueHashtableFromMatrix(cmr, columnHashtable, listmatrix->columnElements, columnData, numColumns, queue, &queueEnd, false) );
  }

  *pnumReductions = 0;
  if (queueEnd > queueStart || (pviolatorSubmatrix && (numRows + numColumns > 0)))
  {
    /* Create list matrix representation. */
    CMR_CALL( CMRlistmatrixInitializeFromMatrix(cmr, listmatrix, matrix) );

    /* We now start main loop. */
    size_t numRowReductions = 0;
    size_t numColumnReductions = 0;
    CMR_CALL( reduceListMatrix(cmr, listmatrix, rowData, columnData, rowHashtable, columnHashtable, hashVector, queue,
      &queueStart, &queueEnd, queueMemory, reductions, pnumReductions, &numRowReductions, &numColumnReductions) );

    if (stats)
    {
      stats->reduceCount++;
      stats->reduceTime += (clock() - reduceClock) * 1.0 / CLOCKS_PER_SEC;
    }

    /* Extract remaining submatrix. */
    if (preducedSubmatrix)
    {
      CMR_CALL( extractRemainingSubmatrix(cmr, matrix, numRowReductions, numColumnReductions, listmatrix,
        preducedSubmatrix) );
    }

    if ((pviolatorSubmatrix || pseparation) && (*pnumReductions != (matrix->numRows + matrix->numColumns)))
    {
      clock_t nonbinaryClock = 0;
      if (stats)
        nonbinaryClock = clock();

      CMR_CALL( calcBinaryHashFromListMatrix(cmr, listmatrix, rowData, columnData, hashVector) );

      queueStart = 0;
      queueEnd = 0;
      CMR_CALL( CMRlisthashtableFree(cmr, &rowHashtable) );
      CMR_CALL( CMRlisthashtableCreate(cmr, &rowHashtable, nextPower2(numRows), numRows) );
      CMR_CALL( initializeQueueHashtableFromListMatrix(cmr, rowHashtable, listmatrix, rowData, queue, &queueEnd,
        true) );

      CMR_CALL( CMRlisthashtableFree(cmr, &columnHashtable) );
      CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns) );
      CMR_CALL( initializeQueueHashtableFromListMatrix(cmr, columnHashtable, listmatrix, columnData, queue, &queueEnd,
        false) );

      CMR_SUBMAT* violatorSubmatrix;
      CMR_CALL( extractNonbinarySubmatrix(cmr, listmatrix, rowData, columnData, rowHashtable, columnHashtable,
        hashVector, queue, &queueStart, &queueEnd, queueMemory, &violatorSubmatrix) );
      if (violatorSubmatrix && pviolatorSubmatrix)
        *pviolatorSubmatrix = violatorSubmatrix;

      if (stats)
      {
        stats->nonbinaryCount++;
        stats->nonbinaryTime += (clock() - nonbinaryClock) * 1.0 / CLOCKS_PER_SEC;
      }

      if (pviolatorSubmatrix && *pviolatorSubmatrix)
      {
        CMRdbgMsg(2, "Extracted an M_2-submatrix.\n");
      }
      else if (!violatorSubmatrix)
      {
        clock_t wheelClock = 0;
        if (stats)
          wheelClock = clock();

        CMR_CALL( extractWheelSubmatrix(cmr, listmatrix, rowData, columnData, queue, queueMemory,
          numRows - numRowReductions, numColumns - numColumnReductions, pviolatorSubmatrix, pseparation) );

        /* Check whether the rank-1 part also has ternary rank 1. */
        if (pseparation)
        {
          CMRdbgMsg(2, "Checking block of -1/+1s for ternary rank 1.\n");

          CMRdbgMsg(2, "Separation has part 0 (%ldx%ld) and part 1 (%ldx%ld) and ranks %d + %d.\n",
            (*pseparation)->numRows[0], (*pseparation)->numColumns[0], (*pseparation)->numRows[1],
            (*pseparation)->numColumns[1], CMRsepaRankTopRight(*pseparation), CMRsepaRank(*pseparation));

          bool sepaIsTernary;
          CMR_SUBMAT* reducedMatrix = NULL;
          if (preducedSubmatrix)
            reducedMatrix = *preducedSubmatrix;
          else
          {
            CMR_CALL( extractRemainingSubmatrix(cmr, matrix, numRowReductions, numColumnReductions, listmatrix,
              &reducedMatrix) );
          }

          CMR_CALL( CMRsepaCheckTernary(cmr, *pseparation, matrix, reducedMatrix, &sepaIsTernary, pviolatorSubmatrix) );
          
          if (!preducedSubmatrix)
            CMR_CALL( CMRsubmatFree(cmr, &reducedMatrix) );
          
          if (!sepaIsTernary)
            CMR_CALL( CMRsepaFree(cmr, pseparation) );
        }

        if (stats)
        {
          stats->wheelCount++;
          stats->wheelTime += (clock() - wheelClock) * 1.0 / CLOCKS_PER_SEC;
        }
      }
      else
      {
        /* We found a violator but the user doesn't want it. */
        assert(!violatorSubmatrix);
        CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
      }
    }
  }
  else
  {
    CMRdbgMsg(2, "No series/parallel element found.\n");

    if (preducedSubmatrix)
      CMR_CALL( createFullRemainingMatrix(cmr, matrix, preducedSubmatrix) );
  }

  CMR_CALL( CMRlisthashtableFree(cmr, &columnHashtable) );
  CMR_CALL( CMRlisthashtableFree(cmr, &rowHashtable) );

  CMR_CALL( CMRfreeStackArray(cmr, &queue) );
  CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowData) );

  CMR_CALL( CMRlistmatrixFree(cmr, &listmatrix) );

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestTernarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool* pisSeriesParallel,
  CMR_SP_REDUCTION* reductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix,
  CMR_SP_STATISTICS* stats)
{
  assert(cmr);
  assert(matrix);
  assert(reductions || !pnumReductions);
  assert(!reductions || pnumReductions);

  CMR_SP_REDUCTION* localReductions = NULL;
  size_t localNumReductions = 0;
  if (!reductions)
    CMR_CALL( CMRallocStackArray(cmr, &localReductions, matrix->numRows + matrix->numColumns) );

  CMR_CALL( decomposeTernarySeriesParallel(cmr, matrix, reductions ? reductions : localReductions, &localNumReductions,
    preducedSubmatrix, pviolatorSubmatrix, NULL, stats) );

  if (pisSeriesParallel)
    *pisSeriesParallel = (*pnumReductions == matrix->numRows + matrix->numColumns);
  if (reductions)
    *pnumReductions = localNumReductions;
  else
    CMR_CALL( CMRfreeStackArray(cmr, &localReductions) );

  return CMR_OKAY;
}

CMR_ERROR CMRdecomposeTernarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool* pisSeriesParallel,
  CMR_SP_REDUCTION* reductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix,
  CMR_SEPA** pseparation, CMR_SP_STATISTICS* stats)
{
  assert(cmr);
  assert(matrix);
  assert(reductions || !pnumReductions);
  assert(!reductions || pnumReductions);

  CMR_SP_REDUCTION* localReductions = NULL;
  size_t localNumReductions = 0;
  if (!reductions)
    CMR_CALL( CMRallocStackArray(cmr, &localReductions, matrix->numRows + matrix->numColumns) );

  CMR_CALL( decomposeTernarySeriesParallel(cmr, matrix, reductions ? reductions : localReductions,
    &localNumReductions, preducedSubmatrix, pviolatorSubmatrix, pseparation, stats) );

  if (pisSeriesParallel)
    *pisSeriesParallel = (*pnumReductions == matrix->numRows + matrix->numColumns);
  if (reductions)
    *pnumReductions = localNumReductions;
  else
    CMR_CALL( CMRfreeStackArray(cmr, &localReductions) );

  return CMR_OKAY;
}

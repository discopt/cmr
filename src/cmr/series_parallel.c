// #define CMR_DEBUG /* Uncomment to debug this file. */

// TODO: Try to not fill the hashtable in initialScan but just add all elements to the queue.

#include <cmr/series_parallel.h>

#include "env_internal.h"
#include "hashtable.h"
#include "sort.h"

#include <limits.h>
#include <stdint.h>
#include <time.h>

#define RANGE_SIGNED_HASH (LLONG_MAX/2)

typedef enum
{
  REMOVED = 0,
  ZERO = 1,
  BLOCK = 2,
  OTHER = 3
} ElementType;

/**
 * \brief Projects \p value into the range [-RANGE_SIGNED_HASH, +RANGE_SIGNED_HASH] via a modulo computation.
 */

static inline
long long projectSignedHash(long long value)
{
  return ((value + RANGE_SIGNED_HASH - 1) % (2*RANGE_SIGNED_HASH-1)) - (RANGE_SIGNED_HASH-1);
}

/**
 * \brief Nonzero in linked-list representation of matrix.
 */

typedef struct _ListNonzero
{
  struct _ListNonzero* left;  /**< \brief Pointer to previous nonzero in the same row. */
  struct _ListNonzero* right; /**< \brief Pointer to next nonzero in the same row. */
  struct _ListNonzero* above; /**< \brief Pointer to previous nonzero in the same column. */
  struct _ListNonzero* below; /**< \brief Pointer to next nonzero in the same column. */
  size_t row;             /**< \brief Row. */
  size_t column;          /**< \brief Column. */
  char value;             /**< \brief Matrix entry. */
  bool disabled;          /**< \brief Whether this edge is disabled in breadth-first search. */
} ListNonzero;

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
  fprintf(stream, "Search for reductions:     %ld / %f\n", stats->reduceCount, stats->reduceTime);
  fprintf(stream, "Search for wheel matrices: %ld / %f\n", stats->wheelCount, stats->wheelTime);
  fprintf(stream, "Total:                     %ld / %f\n", stats->totalCount, stats->totalTime);

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
  ListNonzero* nonzero  /**< \brief Pointer to nonzero. */
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
 * \brief Compares two nonzeros, first by row (ascending) and then by column (ascending) as a tie-breaker.
 */

static
int compareNonzeros(const void* a, const void* b)
{
  const ListNonzero* nonzero1 = (const ListNonzero*)a;
  const ListNonzero* nonzero2 = (const ListNonzero*)b;
  if (nonzero1->row < nonzero2->row)
    return -1;
  else if (nonzero1->row > nonzero2->row)
    return +1;
  else
    return nonzero1->column - nonzero2->column;
}

/**
 * \brief Algorithm data for each element.
 */

typedef struct
{
  ListNonzero nonzeros;                   /**< \brief Dummy nonzero in that row/column. */
  size_t numNonzeros;                 /**< \brief Number of nonzeros in that row/column. */
  long long hashValue;                /**< \brief Hash value of this element. */
  CMR_LISTHASHTABLE_ENTRY hashEntry;   /**< \brief Entry in row or column hashtable. */
  bool inQueue;                       /**< \brief Whether this element is in the queue. */
  char lastBFS;                       /**< \brief Last breadth-first search that found this node.
                                       **< Is 0 initially, positive for search runs, -1 if marked and -2 for SP-reduced
                                       **< element. */
  size_t distance;                    /**< \brief Distance in breadth-first search. */
  size_t predecessor;                 /**< \brief Index of predecessor element in breadth-first search. */
  bool specialBFS;                    /**< \brief Whether this is a special node in breadith-first search. */
} ElementData;

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
    elementData[e].numNonzeros = 0;
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
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Matrix. */
  ElementData* rowData,       /**< Row element data. */
  ElementData* columnData,    /**< Row element data. */
  long long* hashVector       /**< Hash vector. */
)
{
  assert(cmr);
  assert(matrix);
  assert(rowData);
  assert(columnData);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowStarts[row];
    size_t beyond = row+1 < matrix->numRows ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      char value = matrix->entryValues[e];
      assert(value == 1 || value == -1);

      /* Update row data. */
      rowData[row].numNonzeros++;
      long long newHash = projectSignedHash(rowData[row].hashValue + value * hashVector[column]);
      rowData[row].hashValue  = newHash;

      /* Update column data. */
      columnData[column].numNonzeros++;
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
  CMR* cmr,                   /**< \ref CMR environment. */
  ListNonzero* anchor,        /**< Anchor of list representation. */
  ElementData* rowData,       /**< Row element data. */
  ElementData* columnData,    /**< Row element data. */
  long long* hashVector       /**< Hash vector. */
)
{
  assert(cmr);
  assert(anchor);
  assert(rowData);
  assert(columnData);

  /* Reset hash values. */
  for (ListNonzero* rowHead = anchor->below; rowHead != anchor; rowHead = rowHead->below)
    rowData[rowHead->row].hashValue = 0;
  for (ListNonzero* columnHead = anchor->right; columnHead != anchor; columnHead = columnHead->right)
    columnData[columnHead->column].hashValue = 0;

  for (ListNonzero* rowHead = anchor->below; rowHead != anchor; rowHead = rowHead->below)
  {
    for (ListNonzero* nz = rowHead->right; nz != rowHead; nz = nz->right)
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
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_LISTHASHTABLE* hashtable,  /**< Row or column hashtable. */
  ElementData* data,            /**< Row or column data array. */
  size_t sizeData,              /**< Length of \p data. */
  CMR_ELEMENT* queue,            /**< Queue. */
  size_t* pqueueEnd,            /**< Pointer to end of queue. */
  bool isRow                    /**< Whether we are deadling with rows. */
)
{
  assert(cmr);
  assert(hashtable || sizeData == 0);

  for (size_t i = 0; i < sizeData; ++i)
  {
    CMRdbgMsg(2, "%s %d has %d nonzeros.\n", isRow ? "Row" : "Column", i, data[i].numNonzeros);

    /* Check if it qualifies for addition to the hashtable. */
    if (data[i].numNonzeros > 1)
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
  ListNonzero* anchor,          /**< Anchor of list matrix representation. */
  ElementData* data,            /**< Row or column data array. */
  CMR_ELEMENT* queue,           /**< Queue. */
  size_t* pqueueEnd,            /**< Pointer to end of queue. */
  bool isRow                    /**< Whether we are deadling with rows. */
)
{
  assert(cmr);

  CMRdbgMsg(2, "Initializing queue and hashtable from list representation. Inspecting %s.\n",
    isRow ? "rows" : "columns");

  for (ListNonzero* head = isRow ? anchor->below : anchor->right; head != anchor;
    head = (isRow ? head->below : head->right))
  {
    size_t i = isRow ? head->row : head->column;
    CMRdbgMsg(4, "%s %d has %d nonzeros.\n", isRow ? "Row" : "Column", i, data[i].numNonzeros);

    /* Check if it qualifies for addition to the hashtable. */
    assert(data[i].numNonzeros > 1);

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
 * \brief Initializes the list representation of the matrix.
 */

static
CMR_ERROR initializeListMatrix(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,        /**< Matrix. */
  ListNonzero* anchor,          /**< Anchor of row/column data nonzeros. */
  ListNonzero* nonzeros,        /**< Memory for storing the nonzeros. */
  ElementData* rowData,     /**< Row data. */
  ElementData* columnData,  /**< Column data. */
  bool isSorted             /**< Whether the nonzeros in \p matrix are sorted. */
)
{
  size_t i = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowStarts[row];
    size_t beyond = row+1 < matrix->numRows ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (size_t e = first; e < beyond; ++e)
    {
      nonzeros[i].row = row;
      nonzeros[i].column = matrix->entryColumns[e];
      nonzeros[i].value = matrix->entryValues[e];
      nonzeros[i].disabled = false;
      i++;
    }
  }

  /* If necessary, sort the nonzeros in order to create the linked list. */
  if (!isSorted)
    CMR_CALL( CMRsort(cmr, matrix->numNonzeros, nonzeros, sizeof(ListNonzero), compareNonzeros) );
#if !defined(NDEBUG)
  else
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
    {
      size_t start = matrix->rowStarts[row];
      size_t end = row + 1 == matrix->numRows ? matrix->numNonzeros : matrix->rowStarts[row+1];
      for (size_t i = start + 1; i < end; ++i)
      {
        if (matrix->entryColumns[i-1] > matrix->entryColumns[i])
          isSorted = false;
      }
      if (!isSorted)
      {
        fprintf(stderr, "Row r%ld of input matrix is not sorted!", row+1);
        assert("Matrix was expected to be sorted, but is not!" == 0);
      }
    }
  }
#endif /* !NDEBUG */
  
  /* Initialize linked list for rows. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    rowData[row].nonzeros.left = &rowData[row].nonzeros;
    rowData[row].nonzeros.row = row;
    rowData[row].nonzeros.column = SIZE_MAX;
  }
  if (anchor)
  {
    if (matrix->numRows > 0)
    {
      anchor->below = &rowData[0].nonzeros;
      rowData[0].nonzeros.above = anchor;
      anchor->above = &rowData[matrix->numRows-1].nonzeros;
      rowData[matrix->numRows-1].nonzeros.below = anchor;
      for (size_t row = 1; row < matrix->numRows; ++row)
      {
        rowData[row].nonzeros.above = &rowData[row-1].nonzeros;
        rowData[row-1].nonzeros.below = &rowData[row].nonzeros;
      }
    }
    else
    {
      anchor->below = anchor;
      anchor->above = anchor;
    }
  }

  /* Initialize linked list for columns. */
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    columnData[column].nonzeros.above = &columnData[column].nonzeros;
    columnData[column].nonzeros.column = column;
    columnData[column].nonzeros.row = SIZE_MAX;
  }

  if (anchor)
  {
    if (matrix->numColumns > 0)
    {
      anchor->right = &columnData[0].nonzeros;
      columnData[0].nonzeros.left = anchor;
      anchor->left = &columnData[matrix->numColumns-1].nonzeros;
      columnData[matrix->numColumns-1].nonzeros.right = anchor;
      for (size_t column = 1; column < matrix->numColumns; ++column)
      {
        columnData[column].nonzeros.left = &columnData[column-1].nonzeros;
        columnData[column-1].nonzeros.right = &columnData[column].nonzeros;
      }
    }
    else
    {
      anchor->right = anchor;
      anchor->left = anchor;
    }
  }

  /* Link the lists of nonzeros. */
  for (i = 0; i < matrix->numNonzeros; ++i)
  {
    ListNonzero* nz = &nonzeros[i];

    nz->left = rowData[nonzeros[i].row].nonzeros.left;
    rowData[nz->row].nonzeros.left->right = nz;
    rowData[nz->row].nonzeros.left = nz;

    nz->above = columnData[nonzeros[i].column].nonzeros.above;
    columnData[nz->column].nonzeros.above->below = nz;
    columnData[nz->column].nonzeros.above = nz;
  }

  for (size_t row = 0; row < matrix->numRows; ++row)
    rowData[row].nonzeros.left->right = &rowData[row].nonzeros;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnData[column].nonzeros.above->below = &columnData[column].nonzeros;

  return CMR_OKAY;
}

/**
 * \brief Checks whether the row/column at \p index is a (negated) copy of a row/column stored in the hashtable.
 *
 * Compares the actual vectors by traversing the linked-list representation.
 */

static
size_t findCopy(
  ElementData* data,            /**< Row/column data. */
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
      ListNonzero* nz1 = data[index].nonzeros.right;
      ListNonzero* nz2 = data[collisionIndex].nonzeros.right;
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
      ListNonzero* nz1 = data[index].nonzeros.below;
      ListNonzero* nz2 = data[collisionIndex].nonzeros.below;
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
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_LISTHASHTABLE* hashtable,  /**< Row/column hashtable. */
  long long hashChange,         /**< Modification of the hash value. */
  size_t index,                 /**< Index of row/column. */
  ElementData* indexData,       /**< Row/column data. */
  CMR_ELEMENT* queue,            /**< Queue. */
  size_t* pqueueEnd,            /**< Pointer to end of queue. */
  size_t queueMemory,           /**< Memory allocated for queue. */
  bool isRow                    /**< Whether we are dealing with rows. */
)
{
  assert(cmr);
  assert(indexData);
  assert(hashtable);
  assert(queue);

  indexData->numNonzeros--;
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
  CMR* cmr,                             /**< \ref CMR environment. */
  ListNonzero* anchor,                /**< Anchor nonzero. */
  ElementData* rowData,               /**< Row data. */
  ElementData* columnData,            /**< Column data. */
  CMR_LISTHASHTABLE* rowHashtable,     /**< Row hashtable. */
  CMR_LISTHASHTABLE* columnHashtable,  /**< Column hashtable. */
  long long* entryToHash,             /**< Pre-computed hash values of vector entries. */
  CMR_ELEMENT* queue,                  /**< Queue. */
  size_t* pqueueStart,                /**< Pointer to start of queue. */
  size_t* pqueueEnd,                  /**< Pointer to end of queue. */
  size_t queueMemory,                 /**< Memory allocated for queue. */
  CMR_SP_REDUCTION* operations,        /**< Array for storing the SP-reductions. Must be sufficiently large, e.g., number
                                       **< of rows + number of columns. */
  size_t* pnumOperations,             /**< Pointer for storing the number of SP-reductions. */
  size_t* pnumRowOperations,          /**< Pointer for storing the number of row operations (may be \c NULL). */
  size_t* pnumColumnOperations        /**< Pointer for storing the number of column operations (may be \c NULL). */
)
{
  while (*pqueueEnd > *pqueueStart)
  {
#if defined(CMR_DEBUG)
    CMRdbgMsg(0, "\n");
    if (anchor)
    {
      CMRdbgMsg(4, "Status:\n");
      for (ListNonzero* nz = anchor->below; nz != anchor; nz = nz->below)
      {
        size_t row = nz->row;
        CMRdbgMsg(6, "Row r%d: %d nonzeros, hashed = %s, hash = %ld", row+1, rowData[row].numNonzeros,
          rowData[row].hashEntry == SIZE_MAX ? "NO" : "YES", rowData[row].hashValue);
        if (rowData[row].hashEntry != SIZE_MAX)
        {
          CMRdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", rowData[row].hashEntry,
            CMRlisthashtableHash(rowHashtable, rowData[row].hashEntry),
            CMRlisthashtableValue(rowHashtable, rowData[row].hashEntry));
        }
        CMRdbgMsg(0, "\n");
      }
      for (ListNonzero* nz = anchor->right; nz != anchor; nz = nz->right)
      {
        size_t column = nz->column;
        CMRdbgMsg(6, "Column c%d: %d nonzeros, hashed = %s, hash = %ld", column+1, columnData[column].numNonzeros,
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
      CMRdbgMsg(6, "Queue #%d @ %d: %s %d\n", q, q % queueMemory,
        CMRelementIsRow(e) ? "row" : "column", CMRelementIsRow(e) ? CMRelementToRowIndex(e) : CMRelementToColumnIndex(e));
    }
    CMRdbgMsg(0, "\n");
#endif /* CMR_DEBUG */

    CMR_ELEMENT element = queue[(*pqueueStart) % queueMemory];
    ++(*pqueueStart);

    CMRdbgMsg(2, "Top element is %s %d with %d nonzeros.\n", CMRelementIsRow(element) ? "row" : "column",
      CMRelementIsRow(element) ? CMRelementToRowIndex(element) : CMRelementToColumnIndex(element),
      CMRelementIsRow(element) ? rowData[CMRelementToRowIndex(element)].numNonzeros :
      columnData[CMRelementToColumnIndex(element)].numNonzeros);

    if (CMRelementIsRow(element))
    {
      /* We consider a row. */

      size_t row1 = CMRelementToRowIndex(element);
      if (rowData[row1].numNonzeros > 1)
      {
        rowData[row1].inQueue = false;
        size_t row2 = findCopy(rowData, rowHashtable, row1, true, false);

        if (row2 == SIZE_MAX)
        {
          CMRdbgMsg(6, "No parallel row found. Inserting row %d.\n", row1);
          CMR_CALL( CMRlisthashtableInsert(cmr, rowHashtable, llabs(rowData[row1].hashValue), row1, &rowData[row1].hashEntry) );
        }
        else
        {
          CMRdbgMsg(6, "Row %d is parallel.\n", row2);

          /* We found a parallel row. */
          operations[*pnumOperations].element = CMRrowToElement(row1);
          operations[*pnumOperations].mate = CMRrowToElement(row2);
          (*pnumOperations)++;
          (*pnumRowOperations)++;

          for (ListNonzero* entry = rowData[row1].nonzeros.right; entry != &rowData[row1].nonzeros;
            entry = entry->right)
          {
            CMRdbgMsg(8, "Processing nonzero at column %d.\n", entry->column);

            unlinkNonzero(entry);
            CMR_CALL( processNonzero(cmr, columnHashtable, -entryToHash[entry->row] * entry->value, entry->column,
              &columnData[entry->column], queue, pqueueEnd, queueMemory, false) );
          }
          rowData[row1].numNonzeros = 0;
          rowData[row1].lastBFS = -2;
          if (anchor)
          {
            rowData[row1].nonzeros.above->below = rowData[row1].nonzeros.below;
            rowData[row1].nonzeros.below->above = rowData[row1].nonzeros.above;
          }
          assert(rowData[row1].nonzeros.left == &rowData[row1].nonzeros);
          assert(rowData[row1].nonzeros.right == &rowData[row1].nonzeros);
        }
      }
      else
      {
        /* Zero or unit row vector. */
        CMRdbgMsg(4, "Processing %s row %d.\n", rowData[row1].numNonzeros == 0 ? "zero" : "unit", row1);

        rowData[row1].inQueue = false;
        if (rowData[row1].numNonzeros)
        {
          ListNonzero* entry = rowData[row1].nonzeros.right;
          size_t column = entry->column;

          CMRdbgMsg(4, "Processing unit row %d with 1 in column %d.\n", row1, column);

          unlinkNonzero(entry);
          rowData[row1].numNonzeros--;
          CMR_CALL( processNonzero(cmr, columnHashtable, -entryToHash[entry->row] * entry->value, column,
            &columnData[column], queue, pqueueEnd, queueMemory, false) );
          operations[*pnumOperations].mate = CMRcolumnToElement(column);
        }
        else
          operations[*pnumOperations].mate = 0;
        operations[*pnumOperations].element = element;
        (*pnumOperations)++;
        (*pnumRowOperations)++;
        rowData[row1].lastBFS = -2;
        if (anchor)
        {
          rowData[row1].nonzeros.above->below = rowData[row1].nonzeros.below;
          rowData[row1].nonzeros.below->above = rowData[row1].nonzeros.above;
        }
      }
    }
    else
    {
      /* We consider a column. */

      size_t column1 = CMRelementToColumnIndex(element);
      if (columnData[column1].numNonzeros > 1)
      {
        columnData[column1].inQueue = false;
        size_t column2 = findCopy(columnData, columnHashtable, column1, false, false);

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
          operations[*pnumOperations].element = CMRcolumnToElement(column1);
          operations[*pnumOperations].mate = CMRcolumnToElement(column2);
          (*pnumOperations)++;
          (*pnumColumnOperations)++;
          columnData[column1].lastBFS = -2;

          for (ListNonzero* entry = columnData[column1].nonzeros.below; entry != &columnData[column1].nonzeros;
            entry = entry->below)
          {
            CMRdbgMsg(8, "Processing nonzero at row %d.\n", entry->row);

            unlinkNonzero(entry);
            CMR_CALL( processNonzero(cmr, rowHashtable, -entryToHash[entry->column] * entry->value, entry->row,
              &rowData[entry->row], queue, pqueueEnd, queueMemory, true) );
          }
          columnData[column1].numNonzeros = 0;
          if (anchor)
          {
            columnData[column1].nonzeros.left->right = columnData[column1].nonzeros.right;
            columnData[column1].nonzeros.right->left = columnData[column1].nonzeros.left;
          }
          assert(columnData[column1].nonzeros.above == &columnData[column1].nonzeros);
          assert(columnData[column1].nonzeros.below == &columnData[column1].nonzeros);
        }
      }
      else
      {
        /* Zero or unit column vector. */
        CMRdbgMsg(4, "Processing %s column %d.\n", columnData[column1].numNonzeros == 0 ? "zero" : "unit", column1);

        columnData[column1].inQueue = false;
        if (columnData[column1].numNonzeros)
        {
          ListNonzero* entry = columnData[column1].nonzeros.below;
          size_t row = entry->row;

          CMRdbgMsg(4, "Processing unit column %d with 1 in row %d.\n", column1, row);

          unlinkNonzero(entry);
          columnData[column1].numNonzeros--;
          CMR_CALL( processNonzero(cmr, rowHashtable, -entryToHash[entry->column] * entry->value, row,
            &rowData[row], queue, pqueueEnd, queueMemory, true) );
          operations[*pnumOperations].mate = CMRrowToElement(row);
        }
        else
          operations[*pnumOperations].mate = 0;
        operations[*pnumOperations].element = element;
        (*pnumOperations)++;
        (*pnumColumnOperations)++;
        columnData[column1].lastBFS = -2;
        if (anchor)
        {
          columnData[column1].nonzeros.left->right = columnData[column1].nonzeros.right;
          columnData[column1].nonzeros.right->left = columnData[column1].nonzeros.left;
        }
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
  ListNonzero* anchor,                /**< Anchor nonzero. */
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
  assert(anchor);
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
      for (ListNonzero* nz = anchor->below; nz != anchor; nz = nz->below)
      {
        size_t row = nz->row;
        CMRdbgMsg(6, "Row r%d: %d nonzeros, hashed = %s, hash = %ld", row+1, rowData[row].numNonzeros,
          rowData[row].hashEntry == SIZE_MAX ? "NO" : "YES", rowData[row].hashValue);
        if (rowData[row].hashEntry != SIZE_MAX)
        {
          CMRdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", rowData[row].hashEntry,
            CMRlisthashtableHash(rowHashtable, rowData[row].hashEntry),
            CMRlisthashtableValue(rowHashtable, rowData[row].hashEntry));
        }
        CMRdbgMsg(0, "\n");
      }
      for (ListNonzero* nz = anchor->right; nz != anchor; nz = nz->right)
      {
        size_t column = nz->column;
        CMRdbgMsg(6, "Column c%d: %d nonzeros, hashed = %s, hash = %ld", column+1, columnData[column].numNonzeros,
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
      CMRelementIsRow(element) ? rowData[CMRelementToRowIndex(element)].numNonzeros :
      columnData[CMRelementToColumnIndex(element)].numNonzeros);

    if (CMRelementIsRow(element))
    {
      /* We consider a row. */

      size_t row1 = CMRelementToRowIndex(element);
      assert(rowData[row1].numNonzeros > 1);
      rowData[row1].inQueue = false;
      size_t row2 = findCopy(rowData, rowHashtable, row1, true, true);

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
        ListNonzero* nz1 = rowData[row1].nonzeros.right;
        ListNonzero* nz2 = rowData[row2].nonzeros.right;
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
      assert(columnData[column1].numNonzeros > 1);
      columnData[column1].inQueue = false;
      size_t column2 = findCopy(columnData, columnHashtable, column1, false, true);

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
        ListNonzero* nz1 = columnData[column1].nonzeros.below;
        ListNonzero* nz2 = columnData[column2].nonzeros.below;
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
  CMR* cmr,                   /**< \ref CMR environment. */
  size_t currentBFS,        /**< Number of this execution of breadth-first search. */
  ElementData* rowData,     /**< Row data. */
  ElementData* columnData,  /**< Column data. */
  CMR_ELEMENT* queue,        /**< Queue. */
  size_t queueMemory,       /**< Memory for queue. */
  CMR_ELEMENT* sources,      /**< Array of source nodes. */
  size_t numSources,        /**< Number of source nodes. */
  CMR_ELEMENT* targets,      /**< Array of target nodes. */
  size_t numTargets,        /**< Number of target nodes. */
  size_t* pfoundTarget,     /**< Pointer for storing the index of the target node found. */
  size_t* pnumEdges         /**< Pointer for storing the number of traversed edges. */
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
      for (ListNonzero* nz = rowData[row].nonzeros.right; nz->column != SIZE_MAX; nz = nz->right)
      {
        /* Skip edge if disabled. */
        if (nz->disabled)
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
      for (ListNonzero* nz = columnData[column].nonzeros.below; nz->row != SIZE_MAX; nz = nz->below)
      {
        /* Skip edge if disabled. */
        if (nz->disabled)
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
  ElementData* rowData,           /**< Row element data. */
  ElementData* columnData,        /**< Column element data. */
  CMR_SUBMAT** preducedSubmatrix  /**< Pointer for storing the reduced submatrix. */
)
{
  assert(cmr);
  assert(matrix);
  assert(preducedSubmatrix);

  CMR_CALL( CMRsubmatCreate(cmr, matrix->numRows - numRowReductions, matrix->numColumns - numColumnReductions,
    preducedSubmatrix) );
  CMR_SUBMAT* remainingSubmatrix = *preducedSubmatrix;
  size_t rowSubmatrix = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (rowData[row].numNonzeros > 0)
      remainingSubmatrix->rows[rowSubmatrix++] = row;
  }
  assert(rowSubmatrix + numRowReductions == matrix->numRows);

  size_t columnSubmatrix = 0;
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    if (columnData[column].numNonzeros > 0)
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
  CMR* cmr,                             /**< \ref CMR environment. */
  size_t numRows,                       /**< Number of rows. */
  size_t numColumns,                    /**< Number of columns. */
  ElementData* rowData,                 /**< Row element data. */
  ElementData* columnData,              /**< Row element data. */
  CMR_ELEMENT* queue,                   /**< Queue memory. */
  size_t queueMemory,                   /**< Size of \p queue. */
  ListNonzero* anchor,                  /**< Anchor of list matrix. */
  CMR_SUBMAT** pwheelSubmatrix,         /**< Pointer for storing the wheel submatrix. */
  CMR_ELEMENT* separationRank1Elements, /**< Array for storing elements of the rank-1 part of a 2-separation. If not
                                         **< \c NULL, it must have sufficient capacity. */
  size_t* pnumSeparationRank1Elements   /**< Pointer for storing the number of elements stored in
                                         **< \p separationRank1Elements (may be \c NULL). */
)
{
  CMRdbgMsg(2, "Searching for wheel graph representation submatrix.\n");

  assert(anchor);
  assert(anchor->below != anchor);
  size_t currentBFS = 0;
  for (size_t row = 0; row < numRows; ++row)
    rowData[row].specialBFS = false;
  for (size_t column = 0; column < numColumns; ++column)
    columnData[column].specialBFS = false;
  ListNonzero** nzBlock = NULL; /* Pointers for simultaneously traversing columns of block. */
  CMR_CALL( CMRallocStackArray(cmr, &nzBlock, numColumns) );
  CMR_ELEMENT* sources = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &sources, numRows) );
  CMR_ELEMENT* targets = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &targets, numColumns) );
  size_t numEdges = 0;
  for (size_t row = 0; row < numRows; ++row)
    numEdges += rowData[row].numNonzeros;
  while (true)
  {
    size_t sourceRow = anchor->below->row;
    rowData[sourceRow].lastBFS = currentBFS;
    rowData[sourceRow].predecessor = SIZE_MAX;
    rowData[sourceRow].distance = 0;
    size_t targetColumn = rowData[sourceRow].nonzeros.right->column;
    rowData[sourceRow].nonzeros.right->disabled = true;

    CMRdbgMsg(4, "Searching for a chordless cycle from r%d to c%d.\n", sourceRow+1, targetColumn+1);

    clock_t t = clock();
    sources[0] = CMRrowToElement(sourceRow);
    targets[0] = CMRcolumnToElement(targetColumn);
    size_t foundTarget = SIZE_MAX;
    currentBFS++;
    CMR_CALL( breadthFirstSearch(cmr, currentBFS, rowData, columnData, queue, queueMemory, sources, 1, targets, 1,
      &foundTarget, 0) );
    rowData[sourceRow].nonzeros.right->disabled = false;
    assert(foundTarget == 0);
    size_t length = columnData[targetColumn].distance + 1;
    fprintf(stderr, "Time for chordless cycle search: %f\n", (clock() - t) * 1.0 / CLOCKS_PER_SEC);

    CMRdbgMsg(4, "Length of cycle is %d.\n", length);

    if (length > 4)
    {
      /* We found a long chordless cycle. Traverse backwards along path and collect rows/columns. */
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
      break;
    }

    /* We have found a 2-by-2 matrix with only 1's. */
    size_t row1 = sourceRow;
    size_t row2 = columnData[targetColumn].predecessor;
    CMRdbgMsg(4, "Growing the 2x2 submatrix with 1's is at r%d, r%d, c%d, c%d.\n", row1+1, row2+1,
      targetColumn+1, rowData[row2].predecessor+1);

    t = clock();

    /* Go trough the two nonzeros of the two rows simultaneously. */
    ListNonzero* nz1 = rowData[row1].nonzeros.right;
    ListNonzero* nz2 = rowData[row2].nonzeros.right;
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
        nzBlock[numTargets] = columnData[nz1->column].nonzeros.below;
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
          nzBlock[j]->disabled = true;
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

    fprintf(stderr, "Time for growing block: %f\n", (clock() - t) * 1.0 / CLOCKS_PER_SEC);

    CMRdbgMsg(4, "Identified %d source rows.\n", numSources);
    assert(numSources >= 2);

    t = clock();

    currentBFS++;
    foundTarget = SIZE_MAX;
    size_t numTraversedEdges = 0;
    CMR_CALL( breadthFirstSearch(cmr, currentBFS, rowData, columnData, queue, queueMemory, sources, numSources,
      targets, numTargets, &foundTarget, &numTraversedEdges) );

    fprintf(stderr, "Time for second search: %f\n", (clock() - t) * 1.0 / CLOCKS_PER_SEC);

    if (foundTarget < SIZE_MAX)
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
        for (ListNonzero* nz = rowData[row2].nonzeros.right; nz->column != SIZE_MAX; nz = nz->right)
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
        for (ListNonzero* nz = columnData[column2].nonzeros.below; nz->row != SIZE_MAX; nz = nz->below)
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

    if (separationRank1Elements && pnumSeparationRank1Elements)
    {
      CMRdbgMsg(4, "No path found. Extracting 2-separation.\n");

      *pnumSeparationRank1Elements = 0;
      /* Collect all rows that are reachable. */
      for (size_t row = 0; row < numRows; ++row)
      {
        if (rowData[row].lastBFS == currentBFS)
          separationRank1Elements[(*pnumSeparationRank1Elements)++] = CMRrowToElement(row);
      }
      /* Collect all columns that are not reachable. */
      for (size_t column = 0; column < numColumns; ++column)
      {
        if (columnData[column].lastBFS >= -1 && columnData[column].lastBFS != currentBFS)
          separationRank1Elements[(*pnumSeparationRank1Elements)++] = CMRcolumnToElement(column);
      }

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
        for (ListNonzero* nz = columnData[column].nonzeros.below; nz->row != SIZE_MAX; nz = nz->below)
        {
          if (t > 0 || !rowData[nz->row].specialBFS)
          {
            rowData[nz->row].numNonzeros--;
            unlinkNonzero(nz);
          }
          else
            nz->disabled = false;
        }
        if (t == 0)
          columnData[column].numNonzeros = numSources;
        else
        {
          columnData[column].numNonzeros = 0;
          columnData[column].nonzeros.left->right = columnData[column].nonzeros.right;
          columnData[column].nonzeros.right->left = columnData[column].nonzeros.left;
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
        for (ListNonzero* nz = rowData[row].nonzeros.right; nz->column != SIZE_MAX; nz = nz->right)
        {
          if (s > 0 || !columnData[nz->column].specialBFS)
          {
            columnData[nz->column].numNonzeros--;
            unlinkNonzero(nz);
          }
          else
            nz->disabled = false;
        }
        if (s == 0)
          rowData[row].numNonzeros = numTargets;
        else
        {
          rowData[row].numNonzeros = 0;
          rowData[row].nonzeros.left->right = rowData[row].nonzeros.right;
          rowData[row].nonzeros.right->left = rowData[row].nonzeros.left;
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
  CMR* cmr,                             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,                   /**< Sparse char matrix. */
  bool isSorted,                        /**< Whether the entries of \p matrix are sorted. */
  CMR_SP_REDUCTION* reductions,         /**< Array for storing the SP-reductions. Must have capacity at least number of
                                         **< rows + number of columns. */
  size_t* pnumReductions,               /**< Pointer for storing the number of SP-reductions. */
  CMR_SUBMAT** preducedSubmatrix,       /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix,      /**< Pointer for storing a wheel-submatrix (may be \c NULL). */
  CMR_ELEMENT* separationRank1Elements, /**< Array for storing elements of the rank-1 part of a 2-separation. If not
                                         **< \c NULL, it must have sufficient capacity. */
  size_t* pnumSeparationRank1Elements,  /**< Pointer for storing the number of elements stored in
                                         **< \p separationRank1Elements (may be \c NULL). */
  CMR_SP_STATISTICS* stats              /**< Pointer to statistics (may be \c NULL). */
)
{
  assert(cmr);
  assert(matrix);
  assert(reductions && pnumReductions);
  assert(separationRank1Elements || !pnumSeparationRank1Elements);
  assert(!separationRank1Elements || pnumSeparationRank1Elements);

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

  /* Initialize element data and hash vector. */
  ElementData* rowData = NULL;
  CMR_CALL( createElementData(cmr, &rowData, matrix->numRows) );
  ElementData* columnData = NULL;
  CMR_CALL( createElementData(cmr, &columnData, matrix->numColumns) );
  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector,
    matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns) );

  /* Scan the matrix to initialize the element data. */
  CMR_CALL( calcNonzeroCountHashFromMatrix(cmr, matrix, rowData, columnData, hashVector) );

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
    CMR_CALL( initializeQueueHashtableFromMatrix(cmr, rowHashtable, rowData, numRows, queue, &queueEnd, true) );
  }
  CMR_LISTHASHTABLE* columnHashtable = NULL;
  if (numColumns > 0)
  {
    CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns) );
    CMR_CALL( initializeQueueHashtableFromMatrix(cmr, columnHashtable, columnData, numColumns, queue, &queueEnd, false) );
  }

  *pnumReductions = 0;
  if (queueEnd > queueStart || (pviolatorSubmatrix && (numRows + numColumns > 0)))
  {
    /* Create nonzeros of matrix. */
    ListNonzero* nonzeros = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &nonzeros, matrix->numNonzeros) );

    /* Dummy in linked lists for row/column data. */
    ListNonzero the_anchor;
    the_anchor.row = SIZE_MAX;
    the_anchor.column = SIZE_MAX;
    ListNonzero* anchor = pviolatorSubmatrix ? &the_anchor : NULL;
    CMR_CALL( initializeListMatrix(cmr, matrix, anchor, nonzeros, rowData, columnData, isSorted) );

    /* We now start main loop. */
    size_t numRowReductions = 0;
    size_t numColumnReductions = 0;
    CMR_CALL( reduceListMatrix(cmr, anchor, rowData, columnData, rowHashtable, columnHashtable, hashVector, queue,
      &queueStart, &queueEnd, queueMemory, reductions, pnumReductions, &numRowReductions, &numColumnReductions) );

    if (stats)
    {
      stats->reduceCount++;
      stats->reduceTime += (clock() - reduceClock) * 1.0 / CLOCKS_PER_SEC;
    }

    /* Extract remaining submatrix. */
    if (preducedSubmatrix)
    {
      CMR_CALL( extractRemainingSubmatrix(cmr, matrix, numRowReductions, numColumnReductions, rowData, columnData,
        preducedSubmatrix) );
    }

    if (pviolatorSubmatrix && (*pnumReductions != (matrix->numRows + matrix->numColumns)))
    {
      clock_t wheelClock = 0;
      if (stats)
        wheelClock = clock();

      CMR_CALL( extractWheelSubmatrix(cmr, matrix->numRows, matrix->numColumns, rowData, columnData, queue, queueMemory,
        anchor, pviolatorSubmatrix, separationRank1Elements, pnumSeparationRank1Elements) );

      if (stats)
      {
        stats->wheelCount++;
        stats->wheelTime += (clock() - wheelClock) * 1.0 / CLOCKS_PER_SEC;
      }
    }

    CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
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

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestBinarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool isSorted, bool* pisSeriesParallel,
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

  CMR_CALL( decomposeBinarySeriesParallel(cmr, matrix, isSorted, reductions ? reductions : localReductions,
    &localNumReductions, preducedSubmatrix, pviolatorSubmatrix, NULL, NULL, stats) );

  if (pisSeriesParallel)
    *pisSeriesParallel = (*pnumReductions == matrix->numRows + matrix->numColumns);
  if (reductions)
    *pnumReductions = localNumReductions;
  else
    CMR_CALL( CMRfreeStackArray(cmr, &localReductions) );

  return CMR_OKAY;
}


CMR_ERROR CMRdecomposeBinarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool isSorted, bool* pisSeriesParallel,
  CMR_SP_REDUCTION* reductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix,
  CMR_ELEMENT* separationRank1Elements, size_t* pnumSeparationRank1Elements, CMR_SP_STATISTICS* stats)
{
  assert(cmr);
  assert(matrix);
  assert(reductions || !pnumReductions);
  assert(!reductions || pnumReductions);
  assert(separationRank1Elements || !pnumSeparationRank1Elements);
  assert(!separationRank1Elements || pnumSeparationRank1Elements);

  CMR_SP_REDUCTION* localReductions = NULL;
  size_t localNumReductions = 0;
  if (!reductions)
    CMR_CALL( CMRallocStackArray(cmr, &localReductions, matrix->numRows + matrix->numColumns) );

  CMR_CALL( decomposeBinarySeriesParallel(cmr, matrix, isSorted, reductions ? reductions : localReductions,
    &localNumReductions, preducedSubmatrix, pviolatorSubmatrix, separationRank1Elements, pnumSeparationRank1Elements,
    stats) );

  if (pisSeriesParallel)
    *pisSeriesParallel = (*pnumReductions == matrix->numRows + matrix->numColumns);
  if (reductions)
    *pnumReductions = localNumReductions;
  else
    CMR_CALL( CMRfreeStackArray(cmr, &localReductions) );

  return CMR_OKAY;
}

static
CMR_ERROR decomposeTernarySeriesParallel(
  CMR* cmr,                             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,                   /**< Sparse char matrix. */
  bool isSorted,                        /**< Whether the entries of \p matrix are sorted. */
  CMR_SP_REDUCTION* reductions,         /**< Array for storing the SP-reductions. Must have capacity at least number of
                                         **< rows + number of columns. */
  size_t* pnumReductions,               /**< Pointer for storing the number of SP-reductions. */
  CMR_SUBMAT** preducedSubmatrix,       /**< Pointer for storing the SP-reduced submatrix (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix,      /**< Pointer for storing a wheel-submatrix (may be \c NULL). */
  CMR_SP_STATISTICS* stats              /**< Pointer to statistics (may be \c NULL). */
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

  /* Initialize element data and hash vector. */
  ElementData* rowData = NULL;
  CMR_CALL( createElementData(cmr, &rowData, matrix->numRows) );
  ElementData* columnData = NULL;
  CMR_CALL( createElementData(cmr, &columnData, matrix->numColumns) );
  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector,
    matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns) );

  /* Scan the matrix to initialize the element data. */
  CMR_CALL( calcNonzeroCountHashFromMatrix(cmr, matrix, rowData, columnData, hashVector) );

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
    CMR_CALL( initializeQueueHashtableFromMatrix(cmr, rowHashtable, rowData, numRows, queue, &queueEnd, true) );
  }
  CMR_LISTHASHTABLE* columnHashtable = NULL;
  if (numColumns > 0)
  {
    CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns) );
    CMR_CALL( initializeQueueHashtableFromMatrix(cmr, columnHashtable, columnData, numColumns, queue, &queueEnd, false) );
  }

  *pnumReductions = 0;
  if (queueEnd > queueStart || (pviolatorSubmatrix && (numRows + numColumns > 0)))
  {
    /* Create nonzeros of matrix. */
    ListNonzero* nonzeros = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &nonzeros, matrix->numNonzeros) );

    /* Dummy in linked lists for row/column data. */
    ListNonzero the_anchor;
    the_anchor.row = SIZE_MAX;
    the_anchor.column = SIZE_MAX;
    ListNonzero* anchor = pviolatorSubmatrix ? &the_anchor : NULL;
    CMR_CALL( initializeListMatrix(cmr, matrix, anchor, nonzeros, rowData, columnData, isSorted) );

    /* We now start main loop. */
    size_t numRowReductions = 0;
    size_t numColumnReductions = 0;
    CMR_CALL( reduceListMatrix(cmr, anchor, rowData, columnData, rowHashtable, columnHashtable, hashVector, queue,
      &queueStart, &queueEnd, queueMemory, reductions, pnumReductions, &numRowReductions, &numColumnReductions) );

    if (stats)
    {
      stats->reduceCount++;
      stats->reduceTime += (clock() - reduceClock) * 1.0 / CLOCKS_PER_SEC;
    }

    /* Extract remaining submatrix. */
    if (preducedSubmatrix)
    {
      CMR_CALL( extractRemainingSubmatrix(cmr, matrix, numRowReductions, numColumnReductions, rowData, columnData,
        preducedSubmatrix) );
    }

    if (pviolatorSubmatrix && (*pnumReductions != (matrix->numRows + matrix->numColumns)))
    {
      clock_t nonbinaryClock = 0;
      if (stats)
        nonbinaryClock = clock();

      CMR_CALL( calcBinaryHashFromListMatrix(cmr, anchor, rowData, columnData, hashVector) );

      queueStart = 0;
      queueEnd = 0;
      CMR_CALL( CMRlisthashtableFree(cmr, &rowHashtable) );
      CMR_CALL( CMRlisthashtableCreate(cmr, &rowHashtable, nextPower2(numRows), numRows) );
      CMR_CALL( initializeQueueHashtableFromListMatrix(cmr, rowHashtable, anchor, rowData, queue, &queueEnd, true) );

      CMR_CALL( CMRlisthashtableFree(cmr, &columnHashtable) );
      CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns) );
      CMR_CALL( initializeQueueHashtableFromListMatrix(cmr, columnHashtable, anchor, columnData, queue, &queueEnd, false) );

      CMR_CALL( extractNonbinarySubmatrix(cmr, anchor, rowData, columnData, rowHashtable, columnHashtable, hashVector,
        queue, &queueStart, &queueEnd, queueMemory, pviolatorSubmatrix) );

      if (stats)
      {
        stats->nonbinaryCount++;
        stats->nonbinaryTime += (clock() - nonbinaryClock) * 1.0 / CLOCKS_PER_SEC;
      }

      if (*pviolatorSubmatrix)
      {
        CMRdbgMsg(2, "Extracted an M_2-submatrix.\n");
      }
      else
      {
        clock_t wheelClock = 0;
        if (stats)
          wheelClock = clock();

        CMR_CALL( extractWheelSubmatrix(cmr, matrix->numRows, matrix->numColumns, rowData, columnData, queue,
          queueMemory, anchor, pviolatorSubmatrix, NULL, NULL) );

        if (stats)
        {
          stats->wheelCount++;
          stats->wheelTime += (clock() - wheelClock) * 1.0 / CLOCKS_PER_SEC;
        }
      }
    }

    CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
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

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestTernarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool isSorted, bool* pisSeriesParallel,
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

  CMR_CALL( decomposeTernarySeriesParallel(cmr, matrix, isSorted, reductions ? reductions : localReductions,
    &localNumReductions, preducedSubmatrix, pviolatorSubmatrix, stats) );

  if (pisSeriesParallel)
    *pisSeriesParallel = (*pnumReductions == matrix->numRows + matrix->numColumns);
  if (reductions)
    *pnumReductions = localNumReductions;
  else
    CMR_CALL( CMRfreeStackArray(cmr, &localReductions) );

  return CMR_OKAY;
}

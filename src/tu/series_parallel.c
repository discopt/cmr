//#define TU_DEBUG /* Uncomment to debug this file. */
// #define TU_DEBUG_MATRIX_LIST /* Uncomment to print linked list of matrix in each iteration. */

#include <tu/series_parallel.h>

#include "env_internal.h"
#include "hashtable.h"
#include "sort.h"

#include <limits.h>
#include <stdint.h>

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

typedef struct _Nonzero
{
  struct _Nonzero* left;  /**< \brief Pointer to previous nonzero in the same row. */
  struct _Nonzero* right; /**< \brief Pointer to next nonzero in the same row. */
  struct _Nonzero* above; /**< \brief Pointer to previous nonzero in the same column. */
  struct _Nonzero* below; /**< \brief Pointer to next nonzero in the same column. */
  size_t row;             /**< \brief Row. */
  size_t column;          /**< \brief Column. */
  char value;             /**< \brief Matrix entry. */
} Nonzero;

static char seriesParallelStringBuffer[32]; /**< Static buffer for \ref TUspString. */

char* TUspString(TU_SP operation, char* buffer)
{
  if (!buffer)
    buffer = seriesParallelStringBuffer;

  if (operation.element > 0)
  {
    if (operation.mate > 0)
      sprintf(buffer, "c%d copy of c%d", operation.element, operation.mate);
    else if (operation.mate < 0)
      sprintf(buffer, "c%d unit at r%d", operation.element, -operation.mate);
    else
      sprintf(buffer, "c%d zero", operation.element);
    return buffer;
  }
  else if (operation.element < 0)
  {
    if (operation.mate > 0)
      sprintf(buffer, "r%d unit at c%d", -operation.element, operation.mate);
    else if (operation.mate < 0)
      sprintf(buffer, "r%d copy of r%d", -operation.element, -operation.mate);
    else
      sprintf(buffer, "r%d zero", -operation.element);
    return buffer;
  }
  else
    return "<invalid series-parallel operations>";
}

/**
 * \brief Removes the given nonzero from the linked-list representation.
 */

static inline
void unlinkNonzero(
  Nonzero* nonzero  /**< \brief Pointer to nonzero. */
)
{
  assert(nonzero);

  TUdbgMsg(4, "Removing (%d,%d) from linked list.\n", nonzero->row, nonzero->column);
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
  const Nonzero* nonzero1 = (const Nonzero*)a;
  const Nonzero* nonzero2 = (const Nonzero*)b;
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
  Nonzero nonzeros;                   /**< Dummy nonzero in that row/column. */
  size_t numNonzeros;                 /**< Number of nonzeros in that row/column. */
  long long hashValue;                /**< Hash value of this element. */               
  TU_LISTHASHTABLE_ENTRY hashEntry;   /**< Entry in row or column hashtable. */
  bool inQueue;                       /**< Whether this element is in the queue. */
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
 * \brief Scans the matrix initially in order to add all rows or columns either to the queue or to the hashtable.
 */

static
TU_ERROR initialScan(
  TU* tu,                       /**< \ref TU environment. */
  TU_LISTHASHTABLE* hashtable,  /**< Row or column hashtable. */
  ElementData* data,            /**< Row or column data array. */
  size_t sizeData,              /**< Length of \p data. */
  TU_ELEMENT* queue,            /**< Queue. */
  size_t* pqueueEnd,            /**< Pointer to end of queue. */
  bool isRow                    /**< Whether we are deadling with rows. */
)
{
  assert(tu);
  assert(hashtable);

  for (size_t i = 0; i < sizeData; ++i)
  {
    TUdbgMsg(2, "%s %d has %d nonzeros.\n", isRow ? "Row" : "Column", i, data[i].numNonzeros);

    /* Check if it qualifies for addition to the hashtable. */
    if (data[i].numNonzeros > 1)
    {
      TU_LISTHASHTABLE_ENTRY entry = TUlisthashtableFindFirst(hashtable, llabs(data[i].hashValue));
      TUdbgMsg(2, "Search for hash %ld of %s %d yields entry %d.\n", data[i].hashOrType, isRow ? "row" : "column", i, entry);
      if (entry == SIZE_MAX)
      {
        TU_CALL( TUlisthashtableInsert(tu, hashtable, llabs(data[i].hashValue), i, &data[i].hashEntry) );
        continue;
      }
    }

    /* If it was not added to the hashtable, we add it to the queue. */
    queue[*pqueueEnd] = isRow ? TUrowToElement(i) : TUcolumnToElement(i);
    data[i].hashEntry = SIZE_MAX;
    data[i].inQueue = true;
    (*pqueueEnd)++;
  }

  return TU_OKAY;
}

/**
 * \brief Initializes the list representation of the matrix.
 */

static
TU_ERROR initListMatrix(
  TU* tu,                   /**< \ref TU environment. */
  TU_CHRMAT* matrix,        /**< Matrix. */
  Nonzero* nonzeros,        /**< Memory for storing the nonzeros. */
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
      i++;
    }
  }

  /* If necessary, sort the nonzeros in order to create the linked list. */
  if (!isSorted)
    TU_CALL( TUsort(tu, matrix->numNonzeros, nonzeros, sizeof(Nonzero), compareNonzeros) );

  /* Initialize heads of linked lists. */
  Nonzero head;
  head.below = &rowData[0].nonzeros;
  head.above = &rowData[matrix->numRows-1].nonzeros;
  head.right = &columnData[0].nonzeros;
  head.left = &columnData[matrix->numColumns-1].nonzeros;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    rowData[row].nonzeros.above = row > 0 ? &rowData[row-1].nonzeros : &head;
    rowData[row].nonzeros.below = row+1 < matrix->numRows ? &rowData[row+1].nonzeros : &head;
    rowData[row].nonzeros.left = &rowData[row].nonzeros;
    rowData[row].nonzeros.row = row;
    rowData[row].nonzeros.column = SIZE_MAX;
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    columnData[column].nonzeros.left = column > 0 ? &columnData[column-1].nonzeros : &head;
    columnData[column].nonzeros.right = column+1 < matrix->numColumns ? &columnData[column+1].nonzeros : &head;
    columnData[column].nonzeros.above = &columnData[column].nonzeros;
    columnData[column].nonzeros.column = column;
    columnData[column].nonzeros.row = SIZE_MAX;
  }

  /* Link the lists of nonzeros. */
  for (i = 0; i < matrix->numNonzeros; ++i)
  {
    Nonzero* nz = &nonzeros[i];

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

  return TU_OKAY;
}

/**
 * \brief Checks whether the row/column at \p index is a (negated) copy of a row/column stored in the hashtable.
 * 
 * Compares the actual vectors by traversing the linked-list representation.
 */

static
size_t findCopy(
  ElementData* data,            /**< Row/column data. */
  TU_LISTHASHTABLE* hashtable,  /**< Row/column hashtable. */
  size_t index,                 /**< Index in \p data. */
  bool isRow                    /**< Whether we are dealing with rows. */
)
{
  TU_LISTHASHTABLE_HASH hash = llabs(data[index].hashValue);
  TUdbgMsg(4, "Processing %s %d with a collision (hash value %ld).\n", isRow ? "row" : "column", index,
    data[index].hashOrType);
  for (TU_LISTHASHTABLE_ENTRY entry = TUlisthashtableFindFirst(hashtable, hash);
    entry != SIZE_MAX; entry = TUlisthashtableFindNext(hashtable, hash, entry))          
  {
    size_t collisionIndex = TUlisthashtableValue(hashtable, entry);
    TUdbgMsg(8, "%s %d has the same hash value. Comparing...\n", isRow ? "Row" : "Column", j);
    bool equal = true;
    bool negated = true;
    if (isRow)
    {
      Nonzero* nz1 = data[index].nonzeros.right;
      Nonzero* nz2 = data[collisionIndex].nonzeros.right;
      while (equal || negated)
      {
        if (nz1->column != nz2->column)
        {
          equal = false;
          negated = false;
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
      Nonzero* nz1 = data[index].nonzeros.below;
      Nonzero* nz2 = data[collisionIndex].nonzeros.below;
      while (equal || negated)
      {
        if (nz1->row != nz2->row)
        {
          equal = false;
          negated = false;
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

    if (equal || negated)
      return collisionIndex;
  }

  return SIZE_MAX;
}

/**
 * \brief Processes the deletion of a nonzero from the linked-list representation.
 */

static
TU_ERROR processNonzero(
  TU* tu,                       /**< \ref TU environment. */
  TU_LISTHASHTABLE* hashtable,  /**< Row/column hashtable. */
  long long hashChange,         /**< Modification of the hash value. */
  size_t index,                 /**< Index of row/column. */
  ElementData* indexData,       /**< Row/column data. */
  TU_ELEMENT* queue,            /**< Queue. */
  size_t* pqueueEnd,            /**< Pointer to end of queue. */
  size_t queueMemory,           /**< Memory allocated for queue. */
  bool isRow                    /**< Whether we are dealing with rows. */
)
{
  assert(tu);
  assert(indexData);
  assert(hashtable);
  assert(queue);

  indexData->numNonzeros--; 
  long long newHash = projectSignedHash(indexData->hashValue + hashChange);
  TUdbgMsg(4, "Processing nonzero. Old hash is %ld, change is %ld, new hash is %ld.\n", indexData->hashOrType, hashChange,
    newHash);
  indexData->hashValue = newHash;

  /* Add to queue if necessary. */
  if (!indexData->inQueue)
  {
    queue[*pqueueEnd % queueMemory] = isRow ? TUrowToElement(index) : TUcolumnToElement(index);
    indexData->inQueue = true;
    (*pqueueEnd)++;
  }

  /* Remove from hashtable if necessary. */
  if (indexData->hashEntry != SIZE_MAX)
  {
    TU_CALL( TUlisthashtableRemove(tu, hashtable, indexData->hashEntry) );
    indexData->hashEntry = SIZE_MAX;
  }

  return TU_OKAY;
}

TU_ERROR TUfindSeriesParallel(TU* tu, TU_CHRMAT* matrix, TU_SP* operations, size_t* pnumOperations,
  TU_SUBMAT** premainingSubmatrix, /*TU_SUBMAT** pwheelSubmatrix, */bool isSorted)
{
  assert(tu);
  assert(matrix);
  assert(operations);
  assert(pnumOperations);
  assert(!premainingSubmatrix || !*premainingSubmatrix);
//   assert(!pwheelSubmatrix || !*pwheelSubmatrix);

  TUdbgMsg(0, "Searching for series/parallel elements in a %dx%d matrix with %d nonzeros.\n", matrix->numRows,
    matrix->numColumns, matrix->numNonzeros);

  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;
  ElementData* rowData = NULL;
  ElementData* columnData = NULL;
  TU_CALL( TUallocStackArray(tu, &rowData, numRows) );
  for (size_t row = 0; row < numRows; ++row)
  {
    rowData[row].numNonzeros = 0;
    rowData[row].hashValue = 0;
    rowData[row].hashEntry = SIZE_MAX;
    rowData[row].inQueue = false;
  }
  TU_CALL( TUallocStackArray(tu, &columnData, numColumns) );
  for (size_t column = 0; column < numColumns; ++column)
  {
    columnData[column].numNonzeros = 0;
    columnData[column].hashValue = 0;
    columnData[column].hashEntry = SIZE_MAX;
    columnData[column].inQueue = false;
  }

  /* We prepare the hashing. Every coordinate has its own value. These are added up for all nonzero entries. */
  long long* entryToHash = NULL;
  size_t numEntries = numRows > numColumns ? numRows : numColumns;
  TU_CALL( TUallocStackArray(tu, &entryToHash, numEntries) );
  size_t h = 1;
  for (size_t e = 0; e < numEntries; ++e)
  {
    entryToHash[e] = h;
    TUdbgMsg(2, "Entry %d has hash %ld.\n", e, h);
    h = projectSignedHash(3 * h);
  }

  /* We scan the matrix once to compute the number of nonzeros and the hash of each row and each column. */
  for (size_t row = 0; row < numRows; ++row)
  {
    size_t first = matrix->rowStarts[row];
    size_t beyond = row+1 < numRows ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      char value = matrix->entryValues[e];
      assert(value == 1 || value == -1);

      /* Update row data. */
      rowData[row].numNonzeros++;
      long long newHash = projectSignedHash(rowData[row].hashValue + value * entryToHash[column]);
      rowData[row].hashValue  = newHash;

      /* Update column data. */
      columnData[column].numNonzeros++;
      newHash = projectSignedHash(columnData[column].hashValue + value * entryToHash[row]);
      columnData[column].hashValue = newHash;
    }
  }

  /* Initialize the queue. */
  TU_ELEMENT* queue = NULL;
  size_t queueStart = 0;
  size_t queueEnd = 0;
  size_t queueMemory = numRows + numColumns;
  TU_CALL( TUallocStackArray(tu, &queue, queueMemory) );

  TU_LISTHASHTABLE* rowHashtable = NULL;
  TU_CALL( TUlisthashtableCreate(tu, &rowHashtable, nextPower2(numRows), numRows) );
  TU_CALL( initialScan(tu, rowHashtable, rowData, numRows, queue, &queueEnd, true) );

  TU_LISTHASHTABLE* columnHashtable = NULL;
  TU_CALL( TUlisthashtableCreate(tu, &columnHashtable, nextPower2(numColumns), numColumns) );
  TU_CALL( initialScan(tu, columnHashtable, columnData, numColumns, queue, &queueEnd, false) );

  if (queueEnd > queueStart/* || (pwheelSubmatrix && numEntries)*/)
  {
    /* Create nonzeros of matrix. */
    Nonzero* nonzeros = NULL;
    TU_CALL( TUallocStackArray(tu, &nonzeros, matrix->numNonzeros) );

    TU_CALL( initListMatrix(tu, matrix, nonzeros, rowData, columnData, isSorted) );

    /* We now start main loop. */
    *pnumOperations = 0;
    size_t numRowOperations = 0;
    size_t numColumnOperations = 0;
    while (queueEnd > queueStart)
    {
      assert(queueEnd - queueStart <= numRows + numColumns);

#if defined(TU_DEBUG_MATRIX_LIST)
      TUdbgMsg(0, "Row-wise matrix via list:\n");
      for (size_t row = 0; row < matrix->numRows; ++row)
      {
        for (Nonzero* nz = rowData[row].nonzeros.right; nz != &rowData[row].nonzeros; nz = nz->right)
        {
          TUdbgMsg(2, "Nonzero at (%d,%d); left is (%d,%d), right is (%d,%d), above is (%d,%d), below is (%d,%d)\n",
            nz->row, nz->column, nz->left->row, nz->left->column, nz->right->row, nz->right->column, nz->above->row,
            nz->above->column, nz->below->row, nz->below->column);
        }
      }
#endif /* TU_DEBUG_MATRIX_LIST */

#if defined(TU_DEBUG)
      TUdbgMsg(0, "\n");
      TUdbgMsg(4, "Status:\n");
      for (size_t row = 0; row < matrix->numRows; ++row)
      {
        TUdbgMsg(6, "Row %d: %d nonzeros, hashed = %s, hash = %ld", row, rowData[row].numNonzeros,
          rowData[row].hashEntry == SIZE_MAX ? "NO" : "YES", rowData[row].hashOrType);
        if (rowData[row].hashEntry != SIZE_MAX)
          TUdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", rowData[row].hashEntry,
            TUlisthashtableHash(rowHashtable, rowData[row].hashEntry),
            TUlisthashtableValue(rowHashtable, rowData[row].hashEntry));
        TUdbgMsg(0, "\n");
      }
      for (size_t column = 0; column < matrix->numColumns; ++column)
      {
        TUdbgMsg(6, "Column %d: %d nonzeros, hashed = %s, hash = %ld", column, columnData[column].numNonzeros,
          columnData[column].hashEntry == SIZE_MAX ? "NO" : "YES", columnData[column].hashOrType);
        if (columnData[column].hashEntry != SIZE_MAX)
          TUdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", columnData[column].hashEntry,
            TUlisthashtableHash(columnHashtable, columnData[column].hashEntry),
            TUlisthashtableValue(columnHashtable, columnData[column].hashEntry));
        TUdbgMsg(0, "\n");
      }
      for (size_t q = queueStart; q < queueEnd; ++q)
      {
        TU_ELEMENT e = queue[q % queueMemory];
        TUdbgMsg(6, "Queue #%d @ %d: %s %d\n", q, q % queueMemory,
          TUelementIsRow(e) ? "row" : "column", TUelementIsRow(e) ? TUelementToRowIndex(e) : TUelementToColumnIndex(e));
      }
      TUdbgMsg(0, "\n");
#endif /* TU_DEBUG */

      TU_ELEMENT element = queue[queueStart % queueMemory];
      ++queueStart;
      
      TUdbgMsg(2, "Top element is %s %d with %d nonzeros.\n", TUelementIsRow(element) ? "row" : "column",
        TUelementIsRow(element) ? TUelementToRowIndex(element) : TUelementToColumnIndex(element),
        TUelementIsRow(element) ? rowData[TUelementToRowIndex(element)].numNonzeros :
        columnData[TUelementToColumnIndex(element)].numNonzeros);

      if (TUelementIsRow(element))
      {
        /* We consider a row. */

        size_t row1 = TUelementToRowIndex(element);
        if (rowData[row1].numNonzeros > 1)
        {
          rowData[row1].inQueue = false;
          size_t row2 = findCopy(rowData, rowHashtable, row1, true);
          
          if (row2 == SIZE_MAX)
          {
            TUdbgMsg(6, "No parallel row found. Inserting row %d.\n", row1);
            TU_CALL( TUlisthashtableInsert(tu, rowHashtable, llabs(rowData[row1].hashValue), row1, &rowData[row1].hashEntry) );
          }
          else
          {
            TUdbgMsg(6, "Row %d is parallel.\n", row2);

            /* We found a parallel row. */
            operations[*pnumOperations].element = TUrowToElement(row1);
            operations[*pnumOperations].mate = TUrowToElement(row2);
            (*pnumOperations)++;
            numRowOperations++;

            for (Nonzero* entry = rowData[row1].nonzeros.right; entry != &rowData[row1].nonzeros;
              entry = entry->right)
            {
              TUdbgMsg(8, "Processing nonzero at column %d.\n", entry->column);

              unlinkNonzero(entry);
              TU_CALL( processNonzero(tu, columnHashtable, -entryToHash[entry->row] * entry->value, entry->column,
                &columnData[entry->column], queue, &queueEnd, queueMemory, false) );
            }
            rowData[row1].numNonzeros = 0;
          }
        }
        else
        {
          /* Zero or unit row vector. */
          TUdbgMsg(4, "Processing %s row %d.\n", rowData[row1].numNonzeros == 0 ? "zero" : "unit", row1);

          rowData[row1].inQueue = false;
          if (rowData[row1].numNonzeros)
          {
            Nonzero* entry = rowData[row1].nonzeros.right;
            size_t column = entry->column;

            TUdbgMsg(4, "Processing unit row %d with 1 in column %d.\n", row1, column);

            unlinkNonzero(entry);
            rowData[row1].numNonzeros--;
            TU_CALL( processNonzero(tu, columnHashtable, -entryToHash[entry->row] * entry->value, column,
              &columnData[column], queue, &queueEnd, queueMemory, false) );
            operations[*pnumOperations].mate = TUcolumnToElement(column);
          }
          else
            operations[*pnumOperations].mate = 0;
          operations[*pnumOperations].element = element;
          (*pnumOperations)++;
          numRowOperations++;
        }
      }
      else
      {
        /* We consider a column. */

        size_t column1 = TUelementToColumnIndex(element);
        if (columnData[column1].numNonzeros > 1)
        {
          columnData[column1].inQueue = false;
          size_t column2 = findCopy(columnData, columnHashtable, column1, false);

          if (column2 == SIZE_MAX)
          {
            TUdbgMsg(6, "No parallel column found. Inserting column %d.\n", column1);
            TU_CALL( TUlisthashtableInsert(tu, columnHashtable, llabs(columnData[column1].hashValue), column1,
              &columnData[column1].hashEntry) );
          }
          else
          {
            TUdbgMsg(6, "Column %d is parallel.\n", column2);

            /* We found a parallel column. */
            operations[*pnumOperations].element = TUcolumnToElement(column1);
            operations[*pnumOperations].mate = TUcolumnToElement(column2);
            (*pnumOperations)++;
            numColumnOperations++;

            for (Nonzero* entry = columnData[column1].nonzeros.below; entry != &columnData[column1].nonzeros;
              entry = entry->below)
            {
              TUdbgMsg(8, "Processing nonzero at row %d.\n", entry->row);

              unlinkNonzero(entry);
              TU_CALL( processNonzero(tu, rowHashtable, -entryToHash[entry->column] * entry->value, entry->row,
                &rowData[entry->row], queue, &queueEnd, queueMemory, true) );
            }
            columnData[column1].numNonzeros = 0;
            assert(columnData[column1].nonzeros.above == &columnData[column1].nonzeros);
            assert(columnData[column1].nonzeros.below == &columnData[column1].nonzeros);
          }
        }
        else
        {
          /* Zero or unit column vector. */
          TUdbgMsg(4, "Processing %s column %d.\n", columnData[column1].numNonzeros == 0 ? "zero" : "unit", column1);

          columnData[column1].inQueue = false;
          if (columnData[column1].numNonzeros)
          {
            Nonzero* entry = columnData[column1].nonzeros.below;
            size_t row = entry->row;

            TUdbgMsg(4, "Processing unit column %d with 1 in row %d.\n", column1, row);

            unlinkNonzero(entry);
            columnData[column1].numNonzeros--;
            TU_CALL( processNonzero(tu, rowHashtable, -entryToHash[entry->column] * entry->value, row,
              &rowData[row], queue, &queueEnd, queueMemory, true) );
            operations[*pnumOperations].mate = TUrowToElement(row);
          }
          else
            operations[*pnumOperations].mate = 0;
          operations[*pnumOperations].element = element;
          (*pnumOperations)++;
          numColumnOperations++;
        }
      }
    }

    /* Extract remaining submatrix. */
    if (premainingSubmatrix)
    {
      TU_CALL( TUsubmatCreate(tu, matrix->numRows - numRowOperations, matrix->numColumns - numColumnOperations,
        premainingSubmatrix) );
      TU_SUBMAT* remainingSubmatrix = *premainingSubmatrix;
      size_t rowSubmatrix = 0;
      for (size_t row = 0; row < matrix->numRows; ++row)
      {
        if (rowData[row].numNonzeros > 0)
          remainingSubmatrix->rows[rowSubmatrix++] = row;
      }
      assert(rowSubmatrix + numRowOperations == matrix->numRows);
      size_t columnSubmatrix = 0;
      for (size_t column = 0; column < matrix->numColumns; ++column)
      {
        if (columnData[column].numNonzeros > 0)
          remainingSubmatrix->columns[columnSubmatrix++] = column;
      }
      assert(columnSubmatrix + numColumnOperations == matrix->numColumns);
    }

    /* Search for a wheel representation submatrix. */
//    if (pwheelSubmatrix && *pnumOperations != numEntries)
//    {
//      TUdbgMsg(2, "Searching for wheel graph representation submatrix.\n");
//
//      for (size_t row = 0; row < numRows; ++row)
//        rowData[row].hashOrType = rowData[row].numNonzeros ? ZERO : REMOVED;
//      for (size_t column = 0; column < numColumns; ++column)
//        columnData[column].hashOrType = columnData[column].numNonzeros ? ZERO : REMOVED;
//
//      /* Find first nonzero row and left-most 1 there. */
//      Nonzero* blockTopLeft = NULL;
//      for (size_t row = 0; row < numRows; ++row)
//      {
//        if (rowData[row].hashOrType != REMOVED)
//        {
//          blockTopLeft = rowData[row].nonzeros.right;
//          break;
//        }
//      }
//      assert(blockTopLeft);
//      rowData[blockTopLeft->row].hashOrType = BLOCK;
//
//      /* Find next 1 to the right of top-left. */
//      Nonzero* blockTopNext = blockTopLeft->right;
//      assert(blockTopLeft->column != SIZE_MAX);
//
//      /* Grow all-1's matrix downwards. */
//      Nonzero* first = blockTopLeft->below;
//      Nonzero* second = blockTopNext->below;
//      Nonzero** blockRowEntry = NULL;
//      TU_CALL( TUallocStackArray(tu, &blockRowEntry, matrix->numRows - numRowOperations) );
//      blockRowEntry[0] = blockTopLeft;
//      size_t numBlockRows = 1;
//      size_t numBlockColumns = 0;
//      while (true)
//      {
//        if (first->row == second->row)
//        {
//          rowData[first->row].hashOrType = BLOCK;
//          blockRowEntry[numBlockRows++] = first;
//          first = first->below;
//          if (first->row == SIZE_MAX)
//            break;
//          second = second->below;
//          if (second->row == SIZE_MAX)
//            break;
//        }
//        else if (first->row < second->row)
//        {
//          first = first->below;
//          if (first->row == SIZE_MAX)
//            break;
//        }
//        else
//        {
//          second = second->below;
//          if (second->row == SIZE_MAX)
//            break;
//        }
//      }
//      
//      TUdbgMsg(2, "Block has grown to %dx2.\n", numBlockRows);
//
//      /* Grow all-1's matrix to the right. */
//      if (numBlockRows == 1)
//      {
//        while (blockTopNext->column != SIZE_MAX)
//        {
//          columnData[blockTopNext->column].hashOrType = BLOCK;
//          blockTopNext = blockTopNext->right;
//          ++numBlockColumns;
//        }
//      }
//      else
//      {
//        size_t maximizer = 0; /* Row (of block) */
//        size_t maxColumn = blockTopLeft->column;
//        for (size_t current = 1; true; current = (current + 1) % numBlockRows)
//        {
//          TUdbgMsg(6, "Current = %d, maximizer = %d, column = %d, maxColumn = %d\n",
//            current, maximizer, blockRowEntry[current]->column, maxColumn);
//          if (current == maximizer)
//          {
//            /* We have found a new column. */
//            assert(columnData[blockRowEntry[current]->column].hashOrType == ZERO);
//            columnData[blockRowEntry[current]->column].hashOrType = BLOCK;
//            ++numBlockColumns;
//            blockRowEntry[current] = blockRowEntry[current]->right;
//            maxColumn = blockRowEntry[current]->column;
//            if (maxColumn == SIZE_MAX)
//              break;
//          }
//          else
//          {
//            /* We scan to the right. */
//            size_t column;
//            while ((column = blockRowEntry[current]->column) < maxColumn)
//              blockRowEntry[current] = blockRowEntry[current]->right;
//
//            if (column == SIZE_MAX)
//              break;
//            else if (column > maxColumn)
//            {
//              maximizer = current;
//              maxColumn = column;
//            }
//          }
//        }
//      }
//      TUdbgMsg(2, "Block has grown to %dx%d.\n", numBlockRows, numBlockColumns);
//
//      /* Check columns parallel to block. */
//      for (size_t row = 0; row < numRows; ++row)
//      {
//        if (rowData[row].hashOrType != BLOCK)
//          continue;
//
//        for (Nonzero* entry = rowData[row].nonzeros.right; entry->column != SIZE_MAX; entry = entry->right)
//        {
//          size_t column = entry->column;
//          if (columnData[column].hashOrType == ZERO)
//            columnData[column].hashOrType = OTHER;
//        }
//      }
//
//      /* Check rows parallel to block. */
//      for (size_t column = 0; column < numColumns; ++column)
//      {
//        if (columnData[column].hashOrType != BLOCK)
//          continue;
//
//        for (Nonzero* entry = columnData[column].nonzeros.below; entry->row != SIZE_MAX; entry = entry->below)
//        {
//          size_t row = entry->row;
//          if (rowData[row].hashOrType == ZERO)
//            rowData[row].hashOrType = OTHER;
//        }
//      }
//
//#if defined(TU_DEBUG)
//      for (size_t row = 0; row < numRows; ++row)
//      {
//        if (rowData[row].hashOrType == BLOCK)
//          TUdbgMsg(4, "Row %d belongs to block.\n", row);
//        else if (rowData[row].hashOrType == ZERO)
//          TUdbgMsg(4, "Row %d is zero parallel to block.\n", row);
//        else if (rowData[row].hashOrType == REMOVED)
//          TUdbgMsg(4, "Row %d was removed.\n", row);
//        else
//          TUdbgMsg(4, "Row %d is nonzero parallel to block.\n", row);
//      }
//      for (size_t column = 0; column < numColumns; ++column)
//      {
//        if (columnData[column].hashOrType == BLOCK)
//          TUdbgMsg(4, "Column %d belongs to block.\n", column);
//        else if (columnData[column].hashOrType == ZERO)
//          TUdbgMsg(4, "Column %d is zero parallel to block.\n", column);
//        else if (columnData[column].hashOrType == REMOVED)
//          TUdbgMsg(4, "Column %d was removed.\n", column);
//        else
//          TUdbgMsg(4, "Column %d is nonzero parallel to block.\n", column);
//      }
//#endif /* TU_DEBUG */
//
//      /* We re-use the queue for a breadth-first search. */
//      queueStart = 0;
//      queueEnd = 0;
//      for (size_t row = 0; row < numRows; ++row)
//      {
//        if (rowData[row].hashOrType == OTHER)
//        {
//          queue[queueEnd++] = TUrowToElement(row);
//          rowData[row].inQueue = true;
//          rowData[row].numNonzeros = SIZE_MAX;
//        }
//      }
//      while (queueStart < queueEnd)
//      {
//        TU_ELEMENT element = queue[queueStart];
//        if (TUelementIsRow(element))
//        {
//          size_t row = TUelementToRowIndex(element);
//          for (Nonzero* entry = rowData[row].nonzeros.right; entry->column != SIZE_MAX; entry = entry->right)
//          {
//            size_t column = entry->column;
//            if (columnData[column].hashOrType == BLOCK)
//              continue;
//
//            columnData[column].numNonzeros = row;
//            if (columnData[column].hashOrType == OTHER)
//            {
//              
//            }
//            else if (!columnData[column].inQueue)
//            {
//              queue[queueEnd++] = TUcolumnToElement(column);
//              columnData[column].inQueue = true;
//            }
//          }
//        }
//      }
//
//      /* Cleanup. */
//      TU_CALL( TUfreeStackArray(tu, &blockRowEntry) );
//    }

    TU_CALL( TUfreeStackArray(tu, &nonzeros) );
  }
  else
  {
    TUdbgMsg(2, "No series/parallel element found.\n");

    if (premainingSubmatrix)
    {
      TU_CALL( TUsubmatCreate(tu, matrix->numRows, matrix->numColumns, premainingSubmatrix) );
      TU_SUBMAT* remainingSubmatrix = *premainingSubmatrix;
      for (size_t row = 0; row < matrix->numRows; ++row)
        remainingSubmatrix->rows[row] = row;
      for (size_t column = 0; column < matrix->numColumns; ++column)
        remainingSubmatrix->columns[column] = column;
    }
  }

  TU_CALL( TUlisthashtableFree(tu, &columnHashtable) );
  TU_CALL( TUlisthashtableFree(tu, &rowHashtable) );

  TU_CALL( TUfreeStackArray(tu, &queue) );
  TU_CALL( TUfreeStackArray(tu, &entryToHash) );
  TU_CALL( TUfreeStackArray(tu, &columnData) );
  TU_CALL( TUfreeStackArray(tu, &rowData) );

  return TU_OKAY;
}

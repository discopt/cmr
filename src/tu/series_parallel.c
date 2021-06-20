// #define TU_DEBUG /* Uncomment to debug this file. */
// #define TU_DEBUG_MATRIX_LIST /* Uncomment to print linked list of matrix in each iteration. */

#include <tu/series_parallel.h>

#include "env_internal.h"
#include "hashtable.h"
#include "sort.h"

#include <limits.h>
#include <stdint.h>

#define RANGE_SIGNED_HASH (LLONG_MAX/2)

static inline
long long projectSignedHash(long long value)
{
  return ((value + RANGE_SIGNED_HASH - 1) % (2*RANGE_SIGNED_HASH-1)) - (RANGE_SIGNED_HASH-1);
}

typedef struct _Nonzero
{
  struct _Nonzero* left;
  struct _Nonzero* right;
  struct _Nonzero* above;
  struct _Nonzero* below;
  size_t row;
  size_t column;
  char value;
} Nonzero;

static inline
void unlinkNonzero(Nonzero* nonzero)
{
  assert(nonzero);
  nonzero->above->below = nonzero->below;
  nonzero->below->above = nonzero->above;
  nonzero->left->right = nonzero->right;
  nonzero->right->left = nonzero->left;
}

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

typedef struct
{
  Nonzero nonzeros;
  size_t numNonzeros;
  long long hash;
  TU_LISTHASHTABLE_ENTRY hashEntry;
  bool inQueue;
} ElementData;

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

static
TU_ERROR initialScan(TU* tu, TU_LISTHASHTABLE* hashtable, ElementData* data, size_t sizeData, TU_ELEMENT* queue,
  size_t* pqueueEnd, bool isRow)
{
  assert(tu);
  assert(hashtable);

  /* Check for unit vector. */
  for (size_t i = 0; i < sizeData; ++i)
  {
    TUdbgMsg(2, "%s %d has %d nonzeros.\n", isRow ? "Row" : "Column", i, data[i].numNonzeros);
    if (data[i].numNonzeros <= 1)
    {
      TUdbgMsg(4, "Initial inspection: %s %s %d found.\n", data[i].numNonzeros == 0 ? "zero" : "unit",
        isRow ? "row" : "column", i);
      queue[*pqueueEnd] = isRow ? TUrowToElement(i) : TUcolumnToElement(i);
      data[i].inQueue = true;
      data[i].hashEntry = SIZE_MAX;
      (*pqueueEnd)++;
    }
    else
    {
      TU_LISTHASHTABLE_ENTRY entry = TUlisthashtableFindFirst(hashtable, llabs(data[i].hash));
      TUdbgMsg(2, "Search for hash %ld of %s %d yields entry %d.\n", data[i].hash, isRow ? "row" : "column", i, entry);
      if (entry == SIZE_MAX)
      {
        TU_CALL( TUlisthashtableInsert(tu, hashtable, llabs(data[i].hash), i, &data[i].hashEntry) );
      }
      else
      {
#if defined(TU_DEBUG)
        size_t other = TUlisthashtableValue(hashtable, entry);
        TUdbgMsg(4, "Initial inspection: potentially parallel %s %d found.\n", isRow ? "row" : "column", other);
#endif /* TU_DEBUG*/
        queue[*pqueueEnd] = isRow ? TUrowToElement(i) : TUcolumnToElement(i);
        data[i].hashEntry = SIZE_MAX;
        data[i].inQueue = true;
        (*pqueueEnd)++;
      }
    }
  }

  return TU_OKAY;
}

static
TU_ERROR initListMatrix(TU* tu, TU_CHRMAT* matrix, Nonzero* nonzeros, ElementData* rowData, ElementData* columnData, bool isSorted)
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

static
size_t findParallel(
  ElementData* data,            /**< Array with row/column data. */
  TU_LISTHASHTABLE* hashtable,  /**< Hash table for row/column hashes. */
  size_t index,                 /**< Index in \p data. */
  bool isRow                    /**< Whether we are dealing with rows. */
)
{
  TU_LISTHASHTABLE_HASH hash = llabs(data[index].hash);
  TUdbgMsg(4, "Processing %s %d with a collision (hash value %ld).\n", isRow ? "row" : "column", index,
    data[index].hash);
  for (TU_LISTHASHTABLE_ENTRY entry = TUlisthashtableFindFirst(hashtable, hash);
    entry != SIZE_MAX; entry = TUlisthashtableFindNext(hashtable, hash, entry))          
  {
    size_t j = TUlisthashtableValue(hashtable, entry);
    TUdbgMsg(8, "%s %d has the same hash value. Comparing...\n", isRow ? "Row" : "Column", j);
    bool equal = true;
    bool negated = true;
    if (isRow)
    {
      Nonzero* nz1 = data[index].nonzeros.right;
      Nonzero* nz2 = data[j].nonzeros.right;
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
      Nonzero* nz2 = data[j].nonzeros.below;
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
      return j;
  }

  return SIZE_MAX;
}

static
TU_ERROR processNonzero(TU* tu, TU_LISTHASHTABLE* hashtable, long long hashChange, size_t index, ElementData* indexData,
  TU_ELEMENT* queue, size_t* pqueueEnd, size_t queueMemory,  bool isRow)
{
  assert(tu);
  assert(indexData);
  assert(hashtable);
  assert(queue);

  indexData->numNonzeros--; 
  long long newHash = projectSignedHash(indexData->hash + hashChange);
  TUdbgMsg(4, "processing nonzero. Old hash is %ld, change is %ld, new hash is %ld.\n", indexData->hash, hashChange,
    newHash);
  indexData->hash = newHash;

  /* Check whether we created a zero or unit element. */
  if (indexData->numNonzeros <= 1)
  {
    TUdbgMsg(10, "Found %s %s %d.\n", indexData->numNonzeros == 0 ? "zero" : "unit", isRow ? "row" : "column", index);

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
  }
  else
  {
    if (indexData->hashEntry != SIZE_MAX)
    {
      TU_CALL( TUlisthashtableRemove(tu, hashtable, indexData->hashEntry) );
      indexData->hashEntry = SIZE_MAX;
    }

    if (!indexData->inQueue)
    {
      queue[*pqueueEnd % queueMemory] = isRow ? TUrowToElement(index) : TUcolumnToElement(index);
      indexData->inQueue = true;
      (*pqueueEnd)++;
    }
  }

  return TU_OKAY;
}

TU_ERROR TUfindSeriesParallel(TU* tu, TU_CHRMAT* matrix, TU_SERIES_PARALLEL* operations, size_t* pnumOperations,
  TU_SUBMAT** premainingSubmatrix, TU_SUBMAT** pwheelSubmatrix, bool isSorted)
{
  assert(tu);
  assert(matrix);
  assert(operations);
  assert(pnumOperations);
  assert(!premainingSubmatrix || !*premainingSubmatrix);
  assert(!pwheelSubmatrix || !*pwheelSubmatrix);

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
    rowData[row].hash = 0;
    rowData[row].hashEntry = SIZE_MAX;
    rowData[row].inQueue = false;
  }
  TU_CALL( TUallocStackArray(tu, &columnData, numColumns) );
  for (size_t column = 0; column < numColumns; ++column)
  {
    columnData[column].numNonzeros = 0;
    columnData[column].hash = 0;
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
      long long newHash = projectSignedHash(rowData[row].hash + value * entryToHash[column]);
      rowData[row].hash  = newHash;

      /* Update column data. */
      columnData[column].numNonzeros++;
      newHash = projectSignedHash(columnData[column].hash + value * entryToHash[row]);
      columnData[column].hash = newHash;
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

  if (queueEnd > queueStart || (pwheelSubmatrix && numEntries))
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
          rowData[row].hashEntry == SIZE_MAX ? "NO" : "YES", rowData[row].hash);
        if (rowData[row].hashEntry != SIZE_MAX)
          TUdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", rowData[row].hashEntry,
            TUlisthashtableHash(rowHashtable, rowData[row].hashEntry),
            TUlisthashtableValue(rowHashtable, rowData[row].hashEntry));
        TUdbgMsg(0, "\n");
      }
      for (size_t column = 0; column < matrix->numColumns; ++column)
      {
        TUdbgMsg(6, "Column %d: %d nonzeros, hashed = %s, hash = %ld", column, columnData[column].numNonzeros,
          columnData[column].hashEntry == SIZE_MAX ? "NO" : "YES", columnData[column].hash);
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
          size_t row2 = findParallel(rowData, rowHashtable, row1, true);
          
          if (row2 == SIZE_MAX)
          {
            TUdbgMsg(6, "No parallel row found. Inserting row %d.\n", row1);
            TU_CALL( TUlisthashtableInsert(tu, rowHashtable, llabs(rowData[row1].hash), row1, &rowData[row1].hashEntry) );
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
          size_t column2 = findParallel(columnData, columnHashtable, column1, false);

          if (column2 == SIZE_MAX)
          {
            TUdbgMsg(6, "No parallel column found. Inserting column %d.\n", column1);
            TU_CALL( TUlisthashtableInsert(tu, columnHashtable, llabs(columnData[column1].hash), column1,
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

// #define TU_DEBUG /* Uncomment to debug this file. */
// #define TU_DEBUG_MATRIX_LIST /* Uncomment to print linked list of matrix in each iteration. */

#include <tu/series_parallel.h>

#include "env_internal.h"
#include "hashtable.h"
#include "sort.h"

#include <stdint.h>

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
  TU_LISTHASHTABLE_HASH hash;
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

TU_ERROR TUfindSeriesParallelBinary(TU* tu, TU_CHRMAT* matrix, TU_SERIES_PARALLEL* operations, size_t* pnumOperations,
  bool isSorted)
{
  assert(tu);
  assert(matrix);
  assert(operations);
  assert(pnumOperations);
  assert(TUisBinaryChr(tu, matrix, NULL));

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
  size_t* entryToHash = NULL;
  size_t numEntries = numRows > numColumns ? numRows : numColumns;
  TU_CALL( TUallocStackArray(tu, &entryToHash, numEntries) );
  size_t h = 1;
  for (size_t e = 0; e < numEntries; ++e)
  {
    entryToHash[e] = h;
    h = 3 * h + 1;
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
      rowData[row].numNonzeros++;
      rowData[row].hash += value * entryToHash[column];
      columnData[column].numNonzeros++;
      columnData[column].hash += value * entryToHash[row];
    }
  }

  /* Initialize the queue. */
  TU_ELEMENT* queue = NULL;
  size_t queueStart = 0;
  size_t queueEnd = 0;
  size_t queueMemory = numRows + numColumns;
  TU_CALL( TUallocStackArray(tu, &queue, queueMemory) );

  /* Check for unit row. */
  for (size_t row = 0; row < numRows; ++row)
  {
    assert(rowData[row].numNonzeros > 0);
    if (rowData[row].numNonzeros == 1)
    {
      TUdbgMsg(4, "Initial inspection: unit row %d found.\n", row);
      queue[queueEnd] = TUrowToElement(row);
      queueEnd++;
      rowData[row].inQueue = true;
    }
  }

  /* Check for unit column. */
  for (size_t column = 0; column < numColumns; ++column)
  {
    assert(columnData[column].numNonzeros > 0);
    if (columnData[column].numNonzeros == 1)
    {
      TUdbgMsg(4, "Initial inspection: unit column %d found.\n", column);
      queue[queueEnd] = TUcolumnToElement(column);
      queueEnd++;
      columnData[column].inQueue = true;
    }
  }

  /* Create and fill row hashtable unless a collision is detected. */
  TU_LISTHASHTABLE* rowHashtable = NULL;
  TU_CALL( TUlisthashtableCreate(tu, &rowHashtable, nextPower2(numRows) << 1, numRows) );
  for (size_t row = 0; row < numRows; ++row)
  {
    /* We skip hashing if unit row. */
    if (rowData[row].inQueue)
      continue;

    TU_LISTHASHTABLE_ENTRY entry = TUlisthashtableFindFirst(rowHashtable, rowData[row].hash);
    TUdbgMsg(2, "Search for hash %d of row %d yields entry %d.\n", rowData[row].hash, row, entry);
    if (entry == SIZE_MAX)
    {
      TU_CALL( TUlisthashtableInsert(tu, rowHashtable, rowData[row].hash, row, &rowData[row].hashEntry) );
    }
    else
    {
#if defined(TU_DEBUG)
      size_t otherRow = TUlisthashtableValue(rowHashtable, entry);
#endif /* TU_DEBUG*/
      TUdbgMsg(4, "Initial inspection: potentially parallel row %d found.\n", otherRow);
      queue[queueEnd] = TUrowToElement(row);
      rowData[row].hashEntry = SIZE_MAX;
      rowData[row].inQueue = true;
      ++queueEnd;
    }
  }

  /* Create and fill column hashtable unless a collision is detected. */
  TU_LISTHASHTABLE* columnHashtable = NULL;
  TU_CALL( TUlisthashtableCreate(tu, &columnHashtable, nextPower2(numColumns) << 1, numColumns) );
  for (size_t column = 0; column < numColumns; ++column)
  {
    /* We skip hashing if unit column. */
    if (columnData[column].inQueue)
      continue;

    TU_LISTHASHTABLE_ENTRY entry = TUlisthashtableFindFirst(columnHashtable, columnData[column].hash);
    if (entry == SIZE_MAX)
    {
      TU_CALL( TUlisthashtableInsert(tu, columnHashtable, columnData[column].hash, column,
        &columnData[column].hashEntry) );
    }
    else
    {
#if defined(TU_DEBUG)
      size_t otherColumn = TUlisthashtableValue(columnHashtable, entry);
#endif /* TU_DEBUG*/
      TUdbgMsg(4, "Initial inspection: potentially parallel column %d found.\n", otherColumn);
      queue[queueEnd] = TUcolumnToElement(column);
      columnData[column].hashEntry = SIZE_MAX;
      columnData[column].inQueue = true;
      ++queueEnd;
    }
  }

  if (queueEnd > queueStart)
  {
    /* There is at least one series or parallel element. */
    TUdbgMsg(2, "There is at least one series/parallel element.\n");

    /* Create nonzeros of matrix. */
    Nonzero* nonzeros = NULL;
    TU_CALL( TUallocStackArray(tu, &nonzeros, matrix->numNonzeros) );
    size_t i = 0;
    for (size_t row = 0; row < numRows; ++row)
    {
      size_t first = matrix->rowStarts[row];
      size_t beyond = row+1 < numRows ? matrix->rowStarts[row+1] : matrix->numNonzeros;
      for (size_t e = first; e < beyond; ++e)
      {
        nonzeros[i].row = row;
        nonzeros[i].column = matrix->entryColumns[e];
        nonzeros[i].value = matrix->entryValues[e];
        i++;
      }
    }
    
    if (!isSorted)
    {
      /* If necessary, sort the nonzeros in order to create the linked list. */
      TU_CALL( TUsort(tu, matrix->numNonzeros, nonzeros, sizeof(Nonzero), compareNonzeros) );
    }

    /* Initialize heads of linked lists. */
    Nonzero head;
    head.below = &rowData[0].nonzeros;
    head.above = &rowData[numRows-1].nonzeros;
    head.right = &columnData[0].nonzeros;
    head.left = &columnData[numColumns-1].nonzeros;
    for (size_t row = 0; row < numRows; ++row)
    {
      rowData[row].nonzeros.above = row > 0 ? &rowData[row-1].nonzeros : &head;
      rowData[row].nonzeros.below = row+1 < numRows ? &rowData[row+1].nonzeros : &head;
      rowData[row].nonzeros.left = &rowData[row].nonzeros;
      rowData[row].nonzeros.right = &rowData[row].nonzeros;
      rowData[row].nonzeros.row = row;
      rowData[row].nonzeros.column = SIZE_MAX;
    }
    for (size_t column = 0; column < numColumns; ++column)
    {
      columnData[column].nonzeros.left = column > 0 ? &columnData[column-1].nonzeros : &head;
      columnData[column].nonzeros.right = column+1 < numColumns ? &columnData[column+1].nonzeros : &head;
      columnData[column].nonzeros.above = &columnData[column].nonzeros;
      columnData[column].nonzeros.below = &columnData[column].nonzeros;
      columnData[column].nonzeros.column = column;
      columnData[column].nonzeros.row = SIZE_MAX;
    }

    /* Link lists. */
    for (i = 0; i < matrix->numNonzeros; ++i)
    {
      nonzeros[i].left = rowData[nonzeros[i].row].nonzeros.left;
      nonzeros[i].right = &rowData[nonzeros[i].row].nonzeros;
      rowData[nonzeros[i].row].nonzeros.left->right = &nonzeros[i];
      rowData[nonzeros[i].row].nonzeros.left = &nonzeros[i];

      nonzeros[i].above = columnData[nonzeros[i].column].nonzeros.above;
      nonzeros[i].below = &columnData[nonzeros[i].column].nonzeros;
      columnData[nonzeros[i].column].nonzeros.above->below = &nonzeros[i];
      columnData[nonzeros[i].column].nonzeros.above = &nonzeros[i];
    }

    /* We now start main loop. */
    *pnumOperations = 0;

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
        TUdbgMsg(6, "Row %d: %d nonzeros, hashed = %s, hash = %d", row, rowData[row].numNonzeros,
          rowData[row].hashEntry == SIZE_MAX ? "NO" : "YES", rowData[row].hash);
        if (rowData[row].hashEntry != SIZE_MAX)
          TUdbgMsg(0, ", hashtable entry: %d with hash=%d, value=%d", rowData[row].hashEntry,
            TUlisthashtableHash(rowHashtable, rowData[row].hashEntry),
            TUlisthashtableValue(rowHashtable, rowData[row].hashEntry));
        TUdbgMsg(0, "\n");
      }
      for (size_t column = 0; column < matrix->numColumns; ++column)
      {
        TUdbgMsg(6, "Column %d: %d nonzeros, hashed = %s, hash = %d", column, columnData[column].numNonzeros,
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

      /* Invariant: each row/column is either 0 or in queue or in hash. */
      for (size_t row = 0; row < numRows; ++row)
      {
        bool inHash = rowData[row].hashEntry != SIZE_MAX;
        bool inQueue = rowData[row].inQueue;
        bool isZero = rowData[row].numNonzeros == 0;
        bool isUnit = rowData[row].numNonzeros == 1;
        if (((inHash ? 1 : 0) + (inQueue ? 1 : 0) + (isZero ? 1 : 0) != 1) || (isUnit && inHash))
        {
          TUdbgMsg(0, "Row %d: inHash = %s, inQueue = %s, isZero = %s, isUnit = %s\n", row, inHash ? "yes" : "no",
            inQueue ? "yes" : "no", isZero ? "yes" : "no", isUnit ? "yes" : "no");
          assert(!"Row invariant not satisfied.");
        }
      }
      for (size_t column = 0; column < numColumns; ++column)
      {
        bool inHash = columnData[column].hashEntry != SIZE_MAX;
        bool inQueue = columnData[column].inQueue;
        bool isZero = columnData[column].numNonzeros == 0;
        bool isUnit = columnData[column].numNonzeros == 1;
        if (((inHash ? 1 : 0) + (inQueue ? 1 : 0) + (isZero ? 1 : 0) != 1) || (isUnit && inHash))
        {
          TUdbgMsg(0, "Column %d: inHash = %s, inQueue = %s, isZero = %s, isUnit = %s\n", column, inHash ? "yes" : "no",
            inQueue ? "yes" : "no", isZero ? "yes" : "no", isUnit ? "yes" : "no");
          assert(!"Column invariant not satisfied.");
        }
      }
#endif /* TU_DEBUG */

      TU_ELEMENT element = queue[queueStart % queueMemory];
      ++queueStart;

      if (TUelementIsRow(element))
      {
        /* We consider a row. */

        size_t row1 = TUelementToRowIndex(element);
        if (rowData[row1].numNonzeros > 1)
        {
          rowData[row1].inQueue = false;
          TU_LISTHASHTABLE_HASH hash = rowData[row1].hash;
          TUdbgMsg(4, "Processing row %d with a collision (hash value %llu).\n", row1, hash);
          size_t row2 = SIZE_MAX;
          for (TU_LISTHASHTABLE_ENTRY entry = TUlisthashtableFindFirst(rowHashtable, hash);
            entry != SIZE_MAX; entry = TUlisthashtableFindNext(rowHashtable, hash, entry))          
          {
            size_t row = TUlisthashtableValue(rowHashtable, entry);
            TUdbgMsg(8, "Row %d has the same hash value. Comparing...\n", row);
            Nonzero* nz1 = rowData[row1].nonzeros.right;
            Nonzero* nz2 = rowData[row].nonzeros.right;
            bool equal = true;
            while (true)
            {
              if (nz1->column != nz2->column)
              {
                equal = false;
                break;
              }
              if (nz1->column == SIZE_MAX)
                break;
              nz1 = nz1->right;
              nz2 = nz2->right;
            }

            if (equal)
            {
              row2 = row;
              break;
            }
          }
          
          if (row2 == SIZE_MAX)
          {
            TUdbgMsg(6, "No parallel row found. Inserting row %d.\n", row1);
            TU_CALL( TUlisthashtableInsert(tu, rowHashtable, rowData[row1].hash, row1, &rowData[row1].hashEntry) );
          }
          else
          {
            TUdbgMsg(6, "Row %d is parallel.\n", row2);

            /* We found a parallel row. */
            operations[*pnumOperations].element = TUrowToElement(row1);
            operations[*pnumOperations].mate = TUrowToElement(row2);
            (*pnumOperations)++;

            for (Nonzero* entry = rowData[row1].nonzeros.right; entry != &rowData[row1].nonzeros;
              entry = entry->right)
            {
              TUdbgMsg(8, "Processing nonzero at column %d.\n", entry->column);

              /* Remove entry from linked list. */
              entry->left->right = entry->right;
              entry->right->left = entry->left;
              entry->above->below = entry->below;
              entry->below->above = entry->above;
              columnData[entry->column].numNonzeros--;
              TU_LISTHASHTABLE_HASH oldHash = columnData[entry->column].hash;
              TU_LISTHASHTABLE_HASH newHash = oldHash - entry->value * entryToHash[entry->row];
              columnData[entry->column].hash = newHash;

              /* Check whether we created a unit column. */
              if (columnData[entry->column].numNonzeros == 1)
              {
                TUdbgMsg(10, "Found unit column %d.\n", entry->column);
                if (!columnData[entry->column].inQueue)
                {
                  /* Add queue. */
                  queue[queueEnd % queueMemory] = TUcolumnToElement(entry->column);
                  columnData[entry->column].inQueue = true;
                  ++queueEnd;
                }

                if (columnData[entry->column].hashEntry != SIZE_MAX)
                {
                  /* Remove from hash. */
                  TU_CALL( TUlisthashtableRemove(tu, columnHashtable, columnData[entry->column].hashEntry) );
                  columnData[entry->column].hashEntry = SIZE_MAX;
                }
              }
              else if (columnData[entry->column].hashEntry != SIZE_MAX)
              {
                /* We may have created a parallel column. We update the hash value. */
                TUdbgMsg(10, "Column %d is hashed. Re-inserting it.", entry->column);

                TU_CALL( TUlisthashtableRemove(tu, columnHashtable, columnData[entry->column].hashEntry) );

                TUdbgMsg(0, " for updated hash %d.\n", newHash);

                /* Check if there is an entry for the updated hash. */
                TU_LISTHASHTABLE_ENTRY hashEntry = TUlisthashtableFindFirst(columnHashtable, newHash );
                if (hashEntry == SIZE_MAX)
                {
                  /* No entry there, so we can insert one. */
                  TU_CALL( TUlisthashtableInsert(tu, columnHashtable, newHash, entry->column,
                    &columnData[entry->column].hashEntry) );
                }
                else
                {
                  /* Otherwise, we add the column to the queue. */
                  queue[queueEnd % queueMemory] = TUcolumnToElement(entry->column);
                  ++queueEnd;
                  columnData[entry->column].inQueue = true;
                  columnData[entry->column].hashEntry = SIZE_MAX;
                }
              }
            }
            rowData[row1].numNonzeros = 0;
          }
        }
        else
        {
          /* Unit row vector. */
          assert(rowData[row1].numNonzeros == 1);

          rowData[row1].inQueue = false;
          Nonzero* entry = rowData[row1].nonzeros.right;
          size_t column = entry->column;

          TUdbgMsg(4, "Processing unit row %d with 1 in column %d.\n", row1, column);

          /* Store operation. */
          operations[*pnumOperations].element = element;
          operations[*pnumOperations].mate = TUcolumnToElement(column);
          (*pnumOperations)++;

          /* Remove entry from linked list. */
          entry->left->right = entry->right;
          entry->right->left = entry->left;
          entry->above->below = entry->below;
          entry->below->above = entry->above;
          rowData[row1].numNonzeros--;
          columnData[column].numNonzeros--;
          TU_LISTHASHTABLE_HASH oldHash = columnData[entry->column].hash;
          TU_LISTHASHTABLE_HASH newHash = oldHash - entry->value * entryToHash[entry->row];
          columnData[entry->column].hash = newHash;

          /* Check whether we created a unit column. */
          if (columnData[column].numNonzeros == 1)
          {
            TUdbgMsg(6, "Found unit column %d.\n", column);
            if (!columnData[column].inQueue)
            {
              queue[queueEnd % queueMemory] = TUcolumnToElement(column);
              columnData[column].inQueue = true;
              ++queueEnd;
            }

            if (columnData[column].hashEntry != SIZE_MAX)
            {
              // TODO: This code is repeated...
              TUdbgMsg(6, "Column %d existed and is now removed from hash.\n", column);
              TU_CALL( TUlisthashtableRemove(tu, columnHashtable, columnData[column].hashEntry) );
              columnData[column].hashEntry = SIZE_MAX;
            }
          }
          else if (columnData[column].hashEntry != SIZE_MAX)
          {
            /* We may have created a parallel column. We update the hash value. */
            TUdbgMsg(6, "Column %d is hashed. Removing it.", column);
            TU_CALL( TUlisthashtableRemove(tu, columnHashtable, columnData[column].hashEntry) );

            /* Check if there is an entry for the updated hash. */
            TU_LISTHASHTABLE_ENTRY hashEntry = TUlisthashtableFindFirst(columnHashtable, newHash);
            if (hashEntry == SIZE_MAX)
            {
              TUdbgMsg(1, "Inserting it with new hash %d.\n", newHash);
              /* No entry there, so we can insert one. */
              TU_CALL( TUlisthashtableInsert(tu, columnHashtable, newHash, column, &columnData[column].hashEntry) );
            }
            else
            {
              TUdbgMsg(1, "Appending it to queue.\n");
              /* Otherwise, we add the column to the queue. */
              queue[queueEnd % queueMemory] = TUcolumnToElement(column);
              columnData[column].hashEntry = SIZE_MAX;
              columnData[column].inQueue = true;
              ++queueEnd;
            }
          }
        }
      }
      else
      {
        /* We consider a column. */

        size_t column1 = TUelementToColumnIndex(element);
        if (columnData[column1].numNonzeros > 1)
        {
          columnData[column1].inQueue = false;
          TU_LISTHASHTABLE_HASH hash = columnData[column1].hash;
          TUdbgMsg(4, "Processing column %d with a collision (hash value %llu).\n", column1, hash);
          size_t column2 = SIZE_MAX;
          for (TU_LISTHASHTABLE_ENTRY entry = TUlisthashtableFindFirst(columnHashtable, hash);
            entry != SIZE_MAX; entry = TUlisthashtableFindNext(columnHashtable, hash, entry))          
          {
            size_t column = TUlisthashtableValue(columnHashtable, entry);
            TUdbgMsg(6, "Column %d has the same hash value. Comparing...\n", column);
            Nonzero* nz1 = columnData[column1].nonzeros.below;
            Nonzero* nz2 = columnData[column].nonzeros.below;
            bool equal = true;
            while (true)
            {
              if (nz1->row != nz2->row)
              {
                equal = false;
                break;
              }
              if (nz1->row == SIZE_MAX)
                break;
              nz1 = nz1->below;
              nz2 = nz2->below;
            }

            if (equal)
            {
              column2 = column;
              break;
            }
          }

          if (column2 == SIZE_MAX)
          {
            TUdbgMsg(6, "No parallel column found. Inserting column %d.\n", column1);
            TU_CALL( TUlisthashtableInsert(tu, columnHashtable, columnData[column1].hash, column1,
              &columnData[column1].hashEntry) );
          }
          else
          {
            TUdbgMsg(6, "Column %d is parallel.\n", column2);

            /* We found a parallel column. */
            operations[*pnumOperations].element = TUcolumnToElement(column1);
            operations[*pnumOperations].mate = TUcolumnToElement(column2);
            (*pnumOperations)++;

            for (Nonzero* entry = columnData[column1].nonzeros.below; entry != &columnData[column1].nonzeros;
              entry = entry->below)
            {
              TUdbgMsg(8, "Processing nonzero at row %d.\n", entry->row);

              /* Remove entry from linked list. */
              entry->left->right = entry->right;
              entry->right->left = entry->left;
              entry->above->below = entry->below;
              entry->below->above = entry->above;
              rowData[entry->row].numNonzeros--;
              TU_LISTHASHTABLE_HASH oldHash = rowData[entry->row].hash;
              TU_LISTHASHTABLE_HASH newHash = oldHash - entry->value * entryToHash[entry->column];
              rowData[entry->row].hash = newHash;

              /* Check whether we created a unit column. */
              if (rowData[entry->row].numNonzeros == 1)
              {
                TUdbgMsg(10, "Found unit row %d.\n", entry->row);
                if (!rowData[entry->row].inQueue)
                {
                  /* Add to queue. */
                  queue[queueEnd % queueMemory] = TUrowToElement(entry->row);
                  rowData[entry->row].inQueue = true;
                  ++queueEnd;
                }
                
                if (rowData[entry->row].hashEntry != SIZE_MAX)
                {
                  /* Remove from hash. */
                  TU_CALL( TUlisthashtableRemove(tu, rowHashtable, rowData[entry->row].hashEntry) );
                  rowData[entry->row].hashEntry = SIZE_MAX;
                }
              }
              else if (rowData[entry->row].hashEntry != SIZE_MAX)
              {
                /* We may have created a parallel row. We update the hash value. */
                TUdbgMsg(10, "Row %d is hashed. Removing it.", entry->row);
                TU_CALL( TUlisthashtableRemove(tu, rowHashtable, rowData[entry->row].hashEntry) );

                /* Check if there is an entry for the updated hash. */
                TU_LISTHASHTABLE_ENTRY hashEntry = TUlisthashtableFindFirst(rowHashtable, newHash);
                if (hashEntry == SIZE_MAX)
                {
                  TUdbgMsg(1, "Inserting it with new hash %d.\n", newHash);
                  /* No entry there, so we can insert one. */
                  TU_CALL( TUlisthashtableInsert(tu, rowHashtable, newHash, entry->row,
                    &rowData[entry->row].hashEntry) );
                }
                else
                {
                  TUdbgMsg(1, "Appending it to queue.\n");
                  /* Otherwise, we add the row to the queue. */
                  queue[queueEnd % queueMemory] = TUrowToElement(entry->row);
                  rowData[entry->row].hashEntry = SIZE_MAX;
                  rowData[entry->row].inQueue = true;
                  ++queueEnd;
                }
              }
            }
            columnData[column1].numNonzeros = 0;
          }
        }
        else
        {
          /* Unit column vector. */
          assert(columnData[column1].numNonzeros == 1);

          columnData[column1].inQueue = false;
          Nonzero* entry = columnData[column1].nonzeros.below;
          size_t row = entry->row;

          TUdbgMsg(4, "Processing unit column %d with 1 in row %d.\n", column1, row);

          /* Store operation. */
          operations[*pnumOperations].element = element;
          operations[*pnumOperations].mate = TUrowToElement(row);
          (*pnumOperations)++;

          /* Remove entry from linked list. */
          entry->left->right = entry->right;
          entry->right->left = entry->left;
          entry->above->below = entry->below;
          entry->below->above = entry->above;
  
          rowData[row].numNonzeros--;
          columnData[column1].numNonzeros--;
          TU_LISTHASHTABLE_HASH oldHash = rowData[entry->row].hash;
          TU_LISTHASHTABLE_HASH newHash = oldHash - entry->value * entryToHash[entry->column];
          rowData[entry->row].hash = newHash;

          /* Remove from column list. */
          Nonzero* left = columnData[column1].nonzeros.left;
          Nonzero* right = columnData[column1].nonzeros.right;
          left->right = right;
          right->left = left;

          /* Check whether we created a unit row. */
          if (rowData[row].numNonzeros == 1)
          {
            TUdbgMsg(6, "Found unit row %d.\n", row);
            if (!rowData[row].inQueue)
            {
              queue[queueEnd % queueMemory] = TUrowToElement(row);
              rowData[row].inQueue = true;
              ++queueEnd;
            }
            
            if (rowData[row].hashEntry != SIZE_MAX)
            {
              // TODO: This code is repeated...
              TUdbgMsg(6, "Row %d existed and is now removed from hash.\n", row);
              TU_CALL( TUlisthashtableRemove(tu, rowHashtable, rowData[row].hashEntry) );
              rowData[row].hashEntry = SIZE_MAX;
            }
          }
          else if (rowData[row].hashEntry != SIZE_MAX)
          {
            /* We may have created a parallel row. We update the hash value. */
            TUdbgMsg(6, "Row %d is hashed. Re-inserting it.\n", row);
            TU_CALL( TUlisthashtableRemove(tu, rowHashtable, rowData[row].hashEntry) );
            
            /* Check if there is an entry for the updated hash. */
            TU_LISTHASHTABLE_ENTRY hashEntry = TUlisthashtableFindFirst(rowHashtable, newHash);
            if (hashEntry == SIZE_MAX)
            {
              /* No entry there, so we can insert one. */
              TU_CALL( TUlisthashtableInsert(tu, rowHashtable, newHash, row, &rowData[row].hashEntry) );
            }
            else
            {
              /* Otherwise, we add the column to the queue. */
              queue[queueEnd % queueMemory] = TUrowToElement(row);
              rowData[row].hashEntry = SIZE_MAX;
              rowData[row].inQueue = true;
              ++queueEnd;
            }
          }
        }
      }
    }

    TU_CALL( TUfreeStackArray(tu, &nonzeros) );
  }
  else
  {
    TUdbgMsg(2, "No series/parallel element found.\n");
  }

  TU_CALL( TUlisthashtableFree(tu, &columnHashtable) );
  TU_CALL( TUlisthashtableFree(tu, &rowHashtable) );

  TU_CALL( TUfreeStackArray(tu, &queue) );
  TU_CALL( TUfreeStackArray(tu, &entryToHash) );
  TU_CALL( TUfreeStackArray(tu, &columnData) );
  TU_CALL( TUfreeStackArray(tu, &rowData) );

  return TU_OKAY;
}

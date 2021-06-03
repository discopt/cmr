#define TU_DEBUG /** Uncomment to debug this file. */

#include <tu/series_parallel.h>

#include "env_internal.h"
#include "sort.h"

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
  size_t hashValue;
  size_t hashBucket;
  bool hashed;
} ElementData;

typedef struct
{
  bool isRow : 1;
  bool isParallel : 1;
  size_t index;
} QueueData;

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
TU_ERROR removeHashEntry(
  long long* hashtable, /**< A hash table for rows or columns. */
  size_t sizeHashTable, /**< The size of the hash table. */
  ElementData* data,    /**< The corresponding row or column data array. */
  size_t index          /**< The row or column whose entry shall be removed. */
)
{
  assert(hashtable);
  assert(data);

  /* Index will always point to the row/column whose bucket will be filled. */
  size_t value = data[index].hashValue;
  size_t bucket = data[index].hashBucket; /* Bucket will always point to the bucket that is to be filled. */
  size_t lastMatchIndex = sizeHashTable;
  for (size_t b = bucket + 1; b != bucket; b = (b + 1) % sizeHashTable)
  {
    size_t i = hashtable[b];
    if (i >= 0)
    {
      if (data[i].hashValue % sizeHashTable == value % sizeHashTable)
        lastMatchIndex = i;
    }
    else if (lastMatchIndex != sizeHashTable)
    {
      /* We move the bucket `b' to `bucket' and continue searching a replacement of that. */
      hashtable[bucket] = lastMatchIndex;
      data[lastMatchIndex].hashBucket = bucket;
      bucket = b;
      b++;
    }
    else
    {
      /* We did not encounter a bucket that may be copied to `bucket'. */
    }
  }
}

TU_ERROR TUfindSeriesParallelChr(TU* tu, TU_CHRMAT* matrix, TU_SERIES_PARALLEL* operations, size_t* pnumOperations)
{
  assert(tu);
  assert(matrix);
  assert(operations);
  assert(pnumOperations);

  TUdbgMsg(0, "Searching for series/parallel elements in a %dx%d matrix.\n", matrix->numRows, matrix->numColumns);
  
  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;
  ElementData* rowData = NULL;
  ElementData* columnData = NULL;
  TU_CALL( TUallocStackArray(tu, &rowData, numRows) );
  for (size_t row = 0; row < numRows; ++row)
  {
    rowData[row].numNonzeros = 0;
    rowData[row].hashValue = 0;
    rowData[row].hashed = false;
  }
  TU_CALL( TUallocStackArray(tu, &columnData, numColumns) );
  for (size_t column = 0; column < numColumns; ++column)
  {
    columnData[column].numNonzeros = 0;
    columnData[column].hashValue = 0;
    columnData[column].hashed = false;
  }
  
  /* We prepare the hashing. Every coordinate has its own value. These are added up for all nonzero entries. */
  size_t* entriesHash = NULL;
  size_t numEntries = numRows > numColumns ? numRows : numColumns;
  TU_CALL( TUallocStackArray(tu, &entriesHash, numEntries) );
  size_t h = 1;
  for (size_t e = 0; e < numEntries; ++e)
  {
    entriesHash[e] = h;
    h *= 3;
  }

  /* We scan the matrix once to compute the number of nonzeros and the hash of each row and each column. */
  bool found = false;
  for (size_t row = 0; row < numRows; ++row)
  {
    size_t first = matrix->rowStarts[row];
    size_t beyond = row+1 < numRows ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      char value = matrix->entryValues[e];
      rowData[row].numNonzeros++;
      rowData[row].hashValue += value * entriesHash[column];
      columnData[column].numNonzeros++;
      columnData[column].hashValue += value * entriesHash[row];
    }
  }

  /* Initialize the queue. */
  QueueData* queue = NULL;
  size_t queueStart = 0;
  size_t queueEnd = 0;
  TU_CALL( TUallocStackArray(tu, &queue, numRows + numColumns) );

  /* Check for unit row. */
  for (size_t row = 0; row < numRows; ++row)
  {
    assert(rowData[row].numNonzeros > 0);
    if (rowData[row].numNonzeros == 1)
    {
      queue[queueEnd].isRow = true;
      queue[queueEnd].isParallel = false;
      queue[queueEnd].index = row;
      queueEnd++;
    }
  }

  /* Check for unit column. */
  for (size_t column = 0; column < numColumns; ++column)
  {
    assert(columnData[column].numNonzeros > 0);
    if (columnData[column].numNonzeros == 1)
    {
      queue[queueEnd].isRow = false;
      queue[queueEnd].isParallel = false;
      queue[queueEnd].index = column;
      queueEnd++;
    }
  }

  /* Create and fill row hashtable unless a collision is detected. */
  long long* rowHashtable = NULL;
  size_t sizeRowHashtable = nextPower2(numRows << 2);
  TU_CALL( TUallocStackArray(tu, &rowHashtable, sizeRowHashtable) );
  for (size_t e = 0; e < sizeRowHashtable; ++e)
    rowHashtable[e] = -1;
  for (size_t row = 0; row < numRows; ++row)
  {
    size_t bucket = rowData[row].hashValue % sizeRowHashtable;
    rowData[row].hashBucket = bucket;
    if (rowHashtable[bucket] >= 0)
    {
      queue[queueEnd].isRow = true;
      queue[queueEnd].isParallel = true;
      queue[queueEnd].index = row;
      ++queueEnd;
      rowData[row].hashed = false;
    }
    else
    {
      rowHashtable[bucket] = row;
      rowData[row].hashed = true;
    }
  }

  /* Create and fill column hashtable unless collision is detected. */
  long long* columnHashtable = NULL;
  size_t sizeColumnHashtable = nextPower2(numColumns << 2);
  TU_CALL( TUallocStackArray(tu, &columnHashtable, sizeColumnHashtable) );
  for (size_t e = 0; e < sizeColumnHashtable; ++e)
    columnHashtable[e] = -1;

  if (!found)
  {
    for (size_t column = 0; column < numColumns; ++column)
    {
      size_t bucket = columnData[column].hashValue % sizeColumnHashtable;
      columnData[column].hashBucket = bucket;
      if (columnHashtable[bucket] >= 0)
      {
        queue[queueEnd].isRow = false;
        queue[queueEnd].isParallel = true;
        queue[queueEnd].index = column;
        ++queueEnd;
        columnData[column].hashed = false;
      }
      else
      {
        columnHashtable[bucket] = column;
        columnData[column].hashed = true;
      }
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

    /* Sort the nonzeros in order to create the linked list. */
    TU_CALL( TUsort(tu, matrix->numNonzeros, nonzeros, sizeof(Nonzero), compareNonzeros) );

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
    }
    for (size_t column = 0; column < numColumns; ++column)
    {
      columnData[column].nonzeros.left = column > 0 ? &columnData[column-1].nonzeros : &head;
      columnData[column].nonzeros.right = column+1 < numColumns ? &columnData[column+1].nonzeros : &head;
      columnData[column].nonzeros.above = &columnData[column].nonzeros;
      columnData[column].nonzeros.below = &columnData[column].nonzeros;
      columnData[column].nonzeros.column = column;
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
      QueueData job = queue[queueStart];
      ++queueStart;

      if (job.isParallel)
      {
        
      }
      else
      {
        if (job.isRow)
        {
          /* We have a unit row. */
          size_t row = job.index;
          Nonzero* entry = rowData[job.index].nonzeros.right;
          size_t column = entry->column;

          TUdbgMsg(4, "Processing unit row %d with 1 in column %d.\n", row, column);

          /* Store operation. */
          operations[*pnumOperations].element = TUrowToElement(row);
          operations[*pnumOperations].mate = TUcolumnToElement(column);
          (*pnumOperations)++;

          /* Remove entry from linked list. */
          entry->left->right = entry->right;
          entry->right->left = entry->left;
          entry->above->below = entry->below;
          entry->below->above = entry->above;
          char value = entry->value;
          rowData[row].numNonzeros--;
          columnData[column].numNonzeros--;

          /* Remove from row list. */
          Nonzero* above = rowData[job.index].nonzeros.above;
          Nonzero* below = rowData[job.index].nonzeros.below;
          above->below = below;
          below->above = above;

          if (columnData[column].numNonzeros == 1)
          {
            /* We may have created a unit column. */
            TUdbgMsg(6, "Found unit column %d.\n", column);
            queue[queueEnd].isRow = false;
            queue[queueEnd].isParallel = false;
            queue[queueEnd].index = column;
            ++queueEnd;
          }
          else if (columnData[column].hashed)
          {
            /* We may have created a parallel column. We update the hash value. */
            TUdbgMsg(6, "Found parallel column %d.\n", column);
            
            /* First remove the hashmap entry. */
            TU_CALL( removeHashEntry(columnHashtable, sizeColumnHashtable, columnData, column) );

            /* Then compute the new hash. */
            columnData[column].hashed = false;
            columnData[column].hashValue -= value * entriesHash[row];

            /* Add it to the queue. */
            queue[queueEnd].isRow = false;
            queue[queueEnd].isParallel = true;
            queue[queueEnd].index = column;
            ++queueEnd;
          }
        }
        else
        {
          
        }
      }
      
      
//       /* Search for unit row vectors. */
//       if (mayHaveRowUnit)
//       {
//         for (Nonzero* nonzero = head.below; nonzero != &head; nonzero = nonzero->below)
//         {
//           TUdbgMsg(4, "Row %d.\n", nonzero->row);
// 
//           size_t row = nonzero->row;
//           assert(rowData[row].numNonzeros > 0);
//           if (rowData[row].numNonzeros == 1)
//           {
//             Nonzero* entry = nonzero->right;
//             size_t column = entry->column;
//             TUdbgMsg(4, "Found unit row %d with 1 in column %d.\n", row, column);
// 
//             /* Store operation. */
//             operations[*pnumOperations].element = TUrowToElement(row);
//             operations[*pnumOperations].mate = TUcolumnToElement(column);
//             (*pnumOperations)++;
// 
//             /* Remove entry from linked list. */
//             entry->left->right = entry->right;
//             entry->right->left = entry->left;
//             entry->above->below = entry->below;
//             entry->below->above = entry->above;
//             rowData[row].numNonzeros--;
//             columnData[column].numNonzeros--;
// 
//             /* Remove from row list. */
//             Nonzero* above = nonzero->above;
//             Nonzero* below = nonzero->below;
//             above->below = below;
//             below->above = above;
//             nonzero = above;
// 
//             /* This may create new series/parallel columns. */
//             mayHaveColumnParallel = true;
//             if (columnData[column].numNonzeros == 1)
//               mayHaveColumnUnit = true;
//           }
//         }
//         mayHaveRowUnit = false;
//       }
    }

    TU_CALL( TUfreeStackArray(tu, &nonzeros) );

    assert(false);
  }
  else
  {
    TUdbgMsg(2, "No series/parallel element found.\n");
  }

  TU_CALL( TUfreeStackArray(tu, &queue) );
  TU_CALL( TUfreeStackArray(tu, &rowHashtable) );
  TU_CALL( TUfreeStackArray(tu, &entriesHash) );
  TU_CALL( TUfreeStackArray(tu, &columnData) );
  TU_CALL( TUfreeStackArray(tu, &rowData) );

  return TU_OKAY;
}

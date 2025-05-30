// #define CMR_DEBUG /* Uncomment to debug this file. */
// #define CMR_DEBUG_MATRICES /* Uncomment to print matrices. */

#include "env_internal.h"
#include "seymour_internal.h"
#include "hashtable.h"
#include "listmatrix.h"

#include <time.h>

typedef struct
{
  long long hashValue;                /**< \brief Hash value of this element. */
  CMR_LISTHASHTABLE_ENTRY hashEntry;  /**< \brief Entry in row or column hashtable. */
} ElementData;

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
#if defined(CMR_DEBUG_REDUCTION)
    CMRdbgMsg(2, "Entry %d has hash %ld.\n", e, h);
#endif /* CMR_DEBUG_REDUCTION */
    h = projectSignedHash(3 * h);
  }

  return CMR_OKAY;
}

/**
 * \brief Scan the matrix to compute the number of nonzeros and the hash of each row and each column.
 */

static
CMR_ERROR calcHashFromMatrix(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,       /**< Matrix. */
  ElementData* rowData,     /**< Row element data. */
  ElementData* columnData,  /**< Column element data. */
  long long* hashVector     /**< Hash vector. */
)
{
  CMR_UNUSED(cmr);

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

      /* Update row data. */
      long long newHash = projectSignedHash(rowData[row].hashValue + hashVector[column]);
      rowData[row].hashValue  = newHash;

      /* Update column data. */
      newHash = projectSignedHash(columnData[column].hashValue + hashVector[row]);
      columnData[column].hashValue = newHash;
    }
  }

  return CMR_OKAY;
}

CMR_ERROR CMRregularitySimpleSearchThreeSeparation(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(queue);

  CMR_SEYMOUR_STATS* stats = task->stats;
  ElementData* rowData = NULL;
  ElementData* columnData = NULL;
  long long* hashVector = NULL;
  CMR_LISTHASHTABLE* rowHashtable = NULL;
  CMR_LISTHASHTABLE* columnHashtable = NULL;

  CMRdbgMsg(6, "Starting search for simple 3-separations for %zux%zu-matrix.\n", task->node->matrix->numRows,
    task->node->matrix->numColumns);
#if defined(CMR_DEBUG_MATRICES)
  CMR_CALL( CMRchrmatPrintDense(cmr, task->node->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG_MATRICES */

  clock_t startTime = clock();
  if (stats)
    stats->simpleThreeSeparationsCount++;

  /* We need the transpose matrix. */
  CMR_SEYMOUR_NODE* node = task->node;
  CMR_CHRMAT* matrix = node->matrix;
  if (!node->transpose)
    CMR_CALL( CMRchrmatTranspose(cmr, matrix, &node->transpose) );
  CMR_CHRMAT* transpose = node->transpose;

  /* Situation A: An entry whose row and column both have only 2 nonzeros:
   *
   * 1 1 0000
   * 1 ? 1100
   * 0 1
   * 0 1
   * 0 0
   * 0 0
   *
   * is clearly a 2-separation whose top-left entry has this property.
   */
  for (size_t row1 = 0; row1 < matrix->numRows; ++row1)
  {
    size_t first = matrix->rowSlice[row1];
    size_t beyond = matrix->rowSlice[row1 + 1];

    if (beyond != first + 2)
      continue;

    size_t* rowColumns = &matrix->entryColumns[first];
    for (size_t c = 0; c < 2; ++c)
    {
      size_t column1 = rowColumns[c];
      size_t numRowColumns = transpose->rowSlice[column1 + 1] - transpose->rowSlice[column1];
      if (numRowColumns != 2)
        continue;

      /* We have found a nonzero at (r,c) such that row r has 2 nonzeros and column c has 2 nonzeros. */
      size_t column2 = rowColumns[1-c];
      size_t row2 = SIZE_MAX;
      for (size_t e = transpose->rowSlice[column1]; e < transpose->rowSlice[column1 + 1]; ++e)
      {
        size_t row = transpose->entryColumns[e];
        if (row != row1)
          row2 = row;
        else
          assert(row == row1);
      }

      CMRdbgMsg(8, "Found a simple 3-separation for rows {r%zu,r%zu} and columns {c%zu,c%zu}.\n", row1+1, row2+1,
        column1+1, column2+1);

      if (stats)
        stats->simpleThreeSeparationsSuccess++;
      CMR_SEPA* separation = NULL;
      CMR_CALL( CMRsepaCreate(cmr, matrix->numRows, matrix->numColumns, &separation) );

      separation->type = CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS;
      for (size_t row = 0; row < matrix->numRows; ++row)
        separation->rowsFlags[row] = (row == row1 || row == row2) ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;
      for (size_t column = 0; column < matrix->numColumns; ++column)
        separation->columnsFlags[column] = (column == column1 || column == column2) ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;

      CMR_SUBMAT* violatorSubmatrix = NULL;
      CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, separation, matrix, transpose, NULL,
        node->isTernary ? &violatorSubmatrix : NULL) );

      if (violatorSubmatrix)
      {
        CMRdbgMsg(8, "-> 2x2 submatrix with bad determinant.\n");

        CMR_CALL( CMRseymourUpdateViolator(cmr, node, violatorSubmatrix) );
        CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      }
      else
      {
        CMR_CALL( CMRregularityDecomposeThreeSum(cmr, task, queue, separation) );
      }

      CMR_CALL( CMRsepaFree(cmr, &separation) );

      goto cleanup;
    }
  }

  /* Situation B: An entry has 2 nonzeros per column and its removal produces duplicate rows. Then this entry, the row
   * with the other entry and the almost-duplicate row form one part.
   *
   * 1 1 1100
   * 1 0 1010
   * 0 1 1100
   * 0 1
   * 0 1
   * 0 0
   * 0 0
   *
   * The transpose version is also considered.
   */

  CMRdbgMsg(2, "Starting search for simple 3-separations based on almost parallel rows/columns.\n");
#if defined(CMR_DEBUG_MATRICES)
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
#endif /* CMR_DEBUG_MATRICES */

  /* Initialize element data and hash vector. */
  CMR_CALL( createElementData(cmr, &rowData, matrix->numRows) );
  CMR_CALL( createElementData(cmr, &columnData, matrix->numColumns) );
  CMR_CALL( createHashVector(cmr, &hashVector,
    matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns) );

  /* Compute hash values for each row and column. */
  CMR_CALL( calcHashFromMatrix(cmr, matrix, rowData, columnData, hashVector) );

  /* Initialize the hashtables. */
  CMR_CALL( CMRlisthashtableCreate(cmr, &rowHashtable, nextPower2(matrix->numRows), matrix->numRows) );
  CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(matrix->numColumns), matrix->numColumns) );

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    CMRdbgMsg(2, "Row r%zu has hash %ld and %zu nonzeros.\n", row+1, rowData[row].hashValue,
      matrix->rowSlice[row+1] - matrix->rowSlice[row]);
    CMR_CALL( CMRlisthashtableInsert(cmr, rowHashtable, llabs(rowData[row].hashValue), row, &rowData[row].hashEntry) );
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    CMRdbgMsg(2, "Column c%zu has hash %ld and %zu nonzeros.\n", column+1, columnData[column].hashValue,
      transpose->rowSlice[column+1] - transpose->rowSlice[column]);
    CMR_CALL( CMRlisthashtableInsert(cmr, columnHashtable, llabs(columnData[column].hashValue), column,
      &columnData[column].hashEntry) );
  }

  /* Go through all candidate nonzeros and see (via hashtable) what happens if one turns it into a zero. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    if (beyond - first == 2)
    {
      for (size_t e = first; e < beyond; ++e)
      {
        size_t column = matrix->entryColumns[e];

        CMRdbgMsg(4, "Considering removal of entry r%zu,c%zu since r%zu has only 2 nonzeros.\n", row+1, column+1,
          row+1);

        long long newHash = projectSignedHash(columnData[column].hashValue - 1 * hashVector[row]);
        CMRdbgMsg(6, "New hash: %ld\n", newHash);
        for (CMR_LISTHASHTABLE_ENTRY entry = CMRlisthashtableFindFirst(columnHashtable, newHash); entry != SIZE_MAX;
          entry = CMRlisthashtableFindNext(columnHashtable, newHash, entry))
        {
          size_t dupColumn = CMRlisthashtableValue(columnHashtable, entry);
          CMRdbgMsg(6, "Hash collision with c%zu.\n", dupColumn+1);

          /* Walk through both columns simultaneously. */
          size_t beyond1 = transpose->rowSlice[column + 1];
          size_t beyond2 = transpose->rowSlice[dupColumn + 1];
          size_t e1 = transpose->rowSlice[column];
          size_t e2 = transpose->rowSlice[dupColumn];
          size_t identicalRow = SIZE_MAX;
          size_t negatedRow = SIZE_MAX;
          bool isDuplicate = true;
          while (e1 < beyond1 || e2 < beyond2)
          {
            size_t r1 = e1 < beyond1 ? transpose->entryColumns[e1] : SIZE_MAX;
            size_t r2 = e2 < beyond2 ? transpose->entryColumns[e2] : SIZE_MAX;
            if (r1 == r2)
            {
              assert(r1 != SIZE_MAX);
              if (transpose->entryValues[e1] == transpose->entryValues[e2])
                identicalRow = r1;
              else
                negatedRow = r1;
              if (identicalRow != SIZE_MAX && negatedRow != SIZE_MAX)
              {
                /* Found a 2x2 violator. */

                CMR_SUBMAT* violatorSubmatrix = NULL;
                CMR_CALL( CMRsubmatCreate2x2(cmr, identicalRow, negatedRow, column, dupColumn, &violatorSubmatrix) );
                CMR_CALL( CMRseymourUpdateViolator(cmr, node, violatorSubmatrix) );
                CMR_CALL( CMRregularityTaskFree(cmr, &task) );

                goto cleanup;
              }
              else
              {
                ++e1;
                ++e2;
              }
            }
            else if (r1 == row)
            {
              assert(r2 != row);
              e1++;
              continue;
            }
            else if (r2 == row)
            {
              assert(r1 != row);
              e2++;
              continue;
            }
            else
            {
              /* Different rows or only one row, but none equal to modified row. */
              isDuplicate = false;
              break;
            }
          }

          if (isDuplicate)
          {
            size_t otherColumn = matrix->entryColumns[e == first ? first+1 : first];

            CMRdbgMsg(8, "Found a simple 3-separation for rows {r%zu} and columns {c%zu,c%zu,c%zu}.\n", row+1, column+1,
              otherColumn+1, dupColumn+1);

            if (stats)
              stats->simpleThreeSeparationsSuccess++;
            CMR_SEPA* separation = NULL;
            CMR_CALL( CMRsepaCreate(cmr, matrix->numRows, matrix->numColumns, &separation) );

            separation->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;
            for (size_t r = 0; r < matrix->numRows; ++r)
              separation->rowsFlags[r] = (r == row) ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;
            for (size_t c = 0; c < matrix->numColumns; ++c)
            {
              separation->columnsFlags[c] = (c == column || c == otherColumn || c == dupColumn)
                ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;
            }

            CMR_SUBMAT* violatorSubmatrix = NULL;
            CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, separation, matrix, transpose, NULL,
              node->isTernary ? &violatorSubmatrix : NULL) );

            if (violatorSubmatrix)
            {
              CMRdbgMsg(8, "-> 2x2 submatrix with bad determinant.\n");
              CMR_CALL( CMRseymourUpdateViolator(cmr, node, violatorSubmatrix) );
              CMR_CALL( CMRregularityTaskFree(cmr, &task) );
            }
            else
              CMR_CALL( CMRregularityDecomposeThreeSum(cmr, task, queue, separation) );

            CMR_CALL( CMRsepaFree(cmr, &separation) );
            goto cleanup;
          }
          else
          {
            CMRdbgMsg(4, "Discarding hash collision.\n");
          }
        }
      }
    }
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    size_t first = transpose->rowSlice[column];
    size_t beyond = transpose->rowSlice[column + 1];
    if (beyond - first == 2)
    {
      for (size_t e = first; e < beyond; ++e)
      {
        size_t row = transpose->entryColumns[e];

        CMRdbgMsg(4, "Considering removal of entry r%zu,c%zu since c%zu has only 2 nonzeros.\n", row+1, column+1,
          column+1);

        long long newHash = projectSignedHash(rowData[row].hashValue - 1 * hashVector[column]);
        CMRdbgMsg(6, "New hash: %ld\n", newHash);
        for (CMR_LISTHASHTABLE_ENTRY entry = CMRlisthashtableFindFirst(rowHashtable, newHash); entry != SIZE_MAX;
          entry = CMRlisthashtableFindNext(rowHashtable, newHash, entry))
        {
          size_t dupRow = CMRlisthashtableValue(rowHashtable, entry);
          CMRdbgMsg(6, "Hash collision with r%zu.\n", dupRow+1);

          /* Walk through both rows simultaneously. */
          size_t beyond1 = matrix->rowSlice[row + 1];
          size_t beyond2 = matrix->rowSlice[dupRow + 1];
          size_t e1 = matrix->rowSlice[row];
          size_t e2 = matrix->rowSlice[dupRow];
          size_t identicalColumn = SIZE_MAX;
          size_t negatedColumn = SIZE_MAX;
          bool isDuplicate = true;
          while (e1 < beyond1 || e2 < beyond2)
          {
            size_t c1 = e1 < beyond1 ? matrix->entryColumns[e1] : SIZE_MAX;
            size_t c2 = e2 < beyond2 ? matrix->entryColumns[e2] : SIZE_MAX;
            if (c1 == c2)
            {
              assert(c1 != SIZE_MAX);
              if (matrix->entryValues[e1] == matrix->entryValues[e2])
                identicalColumn = c1;
              else
                negatedColumn = c1;
              if (identicalColumn != SIZE_MAX && negatedColumn != SIZE_MAX)
              {
                /* Found a 2x2 violator. */

                CMR_SUBMAT* violatorSubmatrix = NULL;
                CMR_CALL( CMRsubmatCreate2x2(cmr, row, dupRow, identicalColumn, negatedColumn, &violatorSubmatrix) );
                CMR_CALL( CMRseymourUpdateViolator(cmr, node, violatorSubmatrix) );
                CMR_CALL( CMRregularityTaskFree(cmr, &task) );

                goto cleanup;
              }
              else
              {
                ++e1;
                ++e2;
              }
            }
            else if (c1 == column)
            {
              assert(c2 != column);
              e1++;
              continue;
            }
            else if (c2 == column)
            {
              assert(c1 != column);
              e2++;
              continue;
            }
            else
            {
              /* Different rows or only one row, but none equal to modified row. */
              isDuplicate = false;
              break;
            }
          }

          if (isDuplicate)
          {
            size_t otherRow = transpose->entryColumns[e == first ? first+1 : first];

            CMRdbgMsg(8, "Found a simple 3-separation for rows {r%zu,r%zu,r%zu} and columns {c%zu}.\n", row+1,
              otherRow+1, dupRow+1, column+1);

            if (stats)
              stats->simpleThreeSeparationsSuccess++;
            CMR_SEPA* separation = NULL;
            CMR_CALL( CMRsepaCreate(cmr, matrix->numRows, matrix->numColumns, &separation) );

            separation->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;
            for (size_t r = 0; r < matrix->numRows; ++r)
              separation->rowsFlags[r] = (r == row || r == otherRow || r == dupRow) ? CMR_SEPA_SECOND : CMR_SEPA_FIRST;
            for (size_t c = 0; c < matrix->numColumns; ++c)
              separation->columnsFlags[c] = (c == column) ? CMR_SEPA_SECOND : CMR_SEPA_FIRST;

            CMR_SUBMAT* violatorSubmatrix = NULL;
            CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, separation, matrix, transpose, NULL,
              node->isTernary ? &violatorSubmatrix : NULL) );

            if (violatorSubmatrix)
            {
              CMRdbgMsg(8, "-> 2x2 submatrix with bad determinant.\n");
              CMR_CALL( CMRseymourUpdateViolator(cmr, node, violatorSubmatrix) );
              CMR_CALL( CMRregularityTaskFree(cmr, &task) );
            }
            else
              CMR_CALL( CMRregularityDecomposeThreeSum(cmr, task, queue, separation) );

            CMR_CALL( CMRsepaFree(cmr, &separation) );
            goto cleanup;
          }
          else
          {
            CMRdbgMsg(4, "Discarding hash collision.\n");
          }
        }
      }
    }
  }

  CMRdbgMsg(2, "No simple 3-separation found.\n");
#if defined(CMR_DEBUG_MATRICES)
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
#endif /* CMR_DEBUG_MATRICES */

  /* We have found nothing, so we note this and re-add the task to the queue. */
  task->node->testedSimpleThreeSeparations = true;
  CMRregularityQueueAdd(queue, task);

cleanup:

  if (rowData)
  {
    CMR_CALL( CMRlisthashtableFree(cmr, &columnHashtable) );
    CMR_CALL( CMRlisthashtableFree(cmr, &rowHashtable) );

    CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );
    CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
    CMR_CALL( CMRfreeStackArray(cmr, &rowData) );
  }

  if (stats)
    stats->simpleThreeSeparationsTime += (clock() - startTime) * 1.0 / CLOCKS_PER_SEC;

  return CMR_OKAY;
}

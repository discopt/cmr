#define CMR_DEBUG /* Uncomment to debug this file. */

#include "regularity_internal.h"

#include "env_internal.h"
#include "matroid_internal.h"

#include <time.h>

/**
 * \brief Element specific data for the enumeration of 3-separations.
 */

typedef struct
{
  short part;       /**< \brief Part this element belongs to; -1 indicates being unassigned. */
  short type[2];    /**< \brief Type of row/column for each part; element lies in span if and only if nonnegative. */
} ElementData;

/**
 * \brief Checks whether the submatrix with rows assigned to \p part and columns assigned to the other has at least
 *        rank 1.
 *
 * If the rank is at least 1, the row and column of a corresponding nonzero are stored in \p rowRepresentative and
 * \p columnRepresentative, respectively.
 */

static
bool findRank1(
  CMR_CHRMAT* matrix,                 /**< Matrix. */
  ElementData* rowData,               /**< Row element data. */
  ElementData* columnData,            /**< Column element data. */
  size_t rowRepresentative[2][2],     /**< Row representatives for each part. */
  size_t columnRepresentative[2][2],  /**< Column representatives for each part. */
  short part                          /**< Part to which the investigated rows belong. */
)
{
  assert(matrix);
  assert(rowData);
  assert(rowRepresentative);
  assert(columnRepresentative);
  assert(part >=0 && part < 2);
  assert(rowRepresentative[part][0] == SIZE_MAX);
  assert(columnRepresentative[1-part][0] == SIZE_MAX);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (rowData[row].part == part)
    {
      size_t first = matrix->rowSlice[row];
      size_t beyond = matrix->rowSlice[row + 1];
      for (size_t entry = first; entry < beyond; ++entry)
      {
        size_t column = matrix->entryColumns[entry];
        if (columnData[column].part == 1-part)
        {
          rowRepresentative[part][0] = row;
          columnRepresentative[1-part][0] = column;
          return true;
        }
      }
    }
  }

  return false;
}

/**
 * \brief Checks whether the submatrix with rows assigned to \p part and columns assigned to the other has at least
 *        rank 2.
 *
 * Assumes that \p rowRepresentative already contains the first nonzero row in that matrix.
 * If the rank is at least 2, \p rowRepresentative and \p columnRepresentative are extended as to indicate a rank-2
 * submatrix.
 */

static
bool findRank2(
  CMR_CHRMAT* matrix,                 /**< Matrix. */
  ElementData* rowData,               /**< Row element data. */
  ElementData* columnData,            /**< Column element data. */
  size_t rowRepresentative[2][2],     /**< Row representatives for each part. */
  size_t columnRepresentative[2][2],  /**< Column representatives for each part. */
  short part                          /**< Part to which the investigated rows belong. */
)
{
  assert(matrix);
  assert(rowData);
  assert(rowRepresentative);
  assert(columnRepresentative);
  assert(part >=0 && part < 2);
  assert(rowRepresentative[part][0] < SIZE_MAX);
  assert(rowRepresentative[part][1] == SIZE_MAX);
  assert(columnRepresentative[1-part][0] < SIZE_MAX);
  assert(columnRepresentative[1-part][1] == SIZE_MAX);

  CMRdbgMsg(14, "findRank2(): first nonzero is in row r%ld.\n", rowRepresentative[part][0]+1);

  for (size_t row = rowRepresentative[part][0] + 1; row < matrix->numRows; ++row)
  {
    if (rowData[row].part == part)
    {
      CMRdbgMsg(14, "findRank2(): processing row r%ld.\n", row+1);
      size_t entry = matrix->rowSlice[row];
      size_t beyond = matrix->rowSlice[row + 1];
      size_t column = (entry < beyond) ? matrix->entryColumns[entry] : SIZE_MAX;
      size_t entryRep = matrix->rowSlice[rowRepresentative[part][0]];
      size_t beyondRep = matrix->rowSlice[rowRepresentative[part][0] + 1];
      size_t columnRep = (entryRep < beyondRep) ? matrix->entryColumns[entryRep] : SIZE_MAX;
      bool isZero = true;
      bool equalRep = true;
      while (column != SIZE_MAX || columnRep != SIZE_MAX)
      {
        CMRdbgMsg(14, "findRank2(): current row's column is c%ld, representative row's column is c%ld.\n", column+1,
          columnRep+1);
        if (column < columnRep)
        {
          CMRdbgMsg(14, "findRank2(): current row has a 1, representative row has a 0.\n");
          if (columnData[column].part == 1-part)
          {
            CMRdbgMsg(14, "findRank2(): inside submatrix!\n");
            /* New row has a 1-entry but representative does not. */
            rowRepresentative[part][1] = row;
            columnRepresentative[1-part][1] = column;
            return true;
          }
          else
          {
            /* We're outside of the submatrix. */
            ++entry;
            column = (entry < beyond) ? matrix->entryColumns[entry] : SIZE_MAX;
          }
        }
        else if (columnRep < column)
        {
          CMRdbgMsg(14, "findRank2(): current row has a 0, representative row has a 1.\n");
          if (columnData[columnRep].part == 1-part)
          {
            CMRdbgMsg(14, "findRank2(): inside submatrix!\n");
            if (!isZero)
            {
              /* We had a common nonzero before, so the 0 here yields rank 2. */
              rowRepresentative[part][1] = row;
              columnRepresentative[1-part][1] = columnRep;
              return true;
            }
            else
              equalRep = false;
          }

          /* We're outside of the submatrix or still a zero vector. */
          ++entryRep;
          columnRep = (entryRep < beyondRep) ? matrix->entryColumns[entryRep] : SIZE_MAX;
        }
        else
        {
          CMRdbgMsg(14, "findRank2(): current row has a 1, representative row has a 1.\n");
          if (columnData[column].part == 1-part)
          {
            CMRdbgMsg(14, "findRank2(): inside submatrix!\n");
            if (!equalRep)
            {
              /* We had a 1 in rep with a 0 here before, so the two 1s here yield rank 2. */
              rowRepresentative[part][1] = row;
              columnRepresentative[1-part][1] = columnRep;
              return true;
            }
            else
              isZero = false;
          }

          /* We're outside of the submatrix. */
          ++entry;
          column = (entry < beyond) ? matrix->entryColumns[entry] : SIZE_MAX;
          ++entryRep;
          columnRep = (entryRep < beyondRep) ? matrix->entryColumns[entryRep] : SIZE_MAX;
        }
      }
    }
  }

  return false;
}

/**
 * \brief Checks whether the submatrix with rows assigned to \p part and columns assigned to the other has at least
 *        rank 3.
 *
 * Assumes that \p rowRepresentative already contains two distinct nonzero rows in that submatrix.
 */

static
bool findRank3(
  CMR_CHRMAT* matrix,                 /**< Matrix. */
  ElementData* rowData,               /**< Row element data. */
  ElementData* columnData,            /**< Column element data. */
  size_t rowRepresentative[2][2],     /**< Row representatives for each part. */
  short part                          /**< Part to which the investigated rows belong. */
)
{
  assert(matrix);
  assert(rowData);
  assert(rowRepresentative);
  assert(part >=0 && part < 2);
  assert(rowRepresentative[part][0] < SIZE_MAX);
  assert(rowRepresentative[part][1] < SIZE_MAX);

  for (size_t row = rowRepresentative[part][1] + 1; row < matrix->numRows; ++row)
  {
    if (rowData[row].part == part)
    {
      CMRdbgMsg(14, "findRank3(): inspecting row r%ld.\n", row+1);
      size_t entry[3] = {
        matrix->rowSlice[row],
        matrix->rowSlice[rowRepresentative[part][0]],
        matrix->rowSlice[rowRepresentative[part][1]]
      };
      size_t beyond[3] = {
        matrix->rowSlice[row + 1],
        matrix->rowSlice[rowRepresentative[part][0] + 1],
        matrix->rowSlice[rowRepresentative[part][1] + 1]
      };
      size_t column[3];
      for (short i = 0; i < 3; ++i)
        column[i] = entry[i] < beyond[i] ? matrix->entryColumns[entry[i]] : SIZE_MAX;

      bool compatible[4] = { true, true, true, true };
      while (column[0] != SIZE_MAX || column[1] != SIZE_MAX || column[2] != SIZE_MAX)
      {
        size_t minColumn = (column[0] < column[1]) ? column[0] : column[1];
        CMRdbgMsg(14, "findRank3(): inspecting column c%ld.\n", minColumn+1);
        if (column[2] < minColumn)
          minColumn = column[2];
        bool nonzero[3] = { column[0] == minColumn, column[1] == minColumn, column[2] == minColumn };

        if (columnData[minColumn].part == 1-part)
        {
          CMRdbgMsg(14, "findRank3(): nonzeros are (%d,%d,%d).\n", nonzero[0] ? 1 : 0, nonzero[1] ? 1 : 0,
            nonzero[2] ? 1 : 0);

          if (nonzero[0])
            compatible[0] = false;
          if (nonzero[0] != nonzero[1])
            compatible[1] = false;
          if (nonzero[0] != nonzero[2])
            compatible[2] = false;
          if (((nonzero[0] ? 1 : 0) + (nonzero[1] ? 1 : 0) + (nonzero[2] ? 1 : 0)) % 2)
            compatible[3] = false;

          CMRdbgMsg(14, "findRank3(): new compatibility is %d,%d,%d,%d.\n", compatible[0] ? 1 : 0, compatible[1] ? 1 : 0,
            compatible[2] ? 1 : 0, compatible[3] ? 1 : 0);
        }

        /* Advance entries with a nonzero. */
        for (short i = 0; i < 3; ++i)
        {
          if (nonzero[i])
          {
            entry[i]++;
            column[i] = entry[i] < beyond[i] ? matrix->entryColumns[entry[i]] : SIZE_MAX;
          }
        }
      }
      if (!compatible[0] && !compatible[1] && !compatible[2] && !compatible[3])
        return true;
    }
  }

  return false;
}


/**
 * \brief Checks whether the given \p row can be attached to \p part.
 */

static
void determineType(
  CMR_CHRMAT* matrix,                 /**< Matrix. */
  ElementData* rowData,               /**< Row element data. */
  ElementData* columnData,            /**< Column element data. */
  size_t rowRepresentative[2][2],     /**< Row representatives for each part. */
  size_t row,                         /**< Row index to check. */
  short part,                         /**< Part to which the investigated rows belong. */
  bool isRow                          /**< Whether we're actually dealing with rows. */
)
{
  assert(matrix);
  assert(rowData);
  assert(rowRepresentative);
  assert(part >=0 && part < 2);

  size_t entry[3] = {
    matrix->rowSlice[row],
    rowRepresentative[part][0] < SIZE_MAX ? matrix->rowSlice[rowRepresentative[part][0]] : SIZE_MAX,
    rowRepresentative[part][1] < SIZE_MAX ? matrix->rowSlice[rowRepresentative[part][1]] : SIZE_MAX
  };
  size_t beyond[3] = {
    matrix->rowSlice[row + 1],
    rowRepresentative[part][0] < SIZE_MAX ? matrix->rowSlice[rowRepresentative[part][0] + 1] : SIZE_MAX,
    rowRepresentative[part][1] < SIZE_MAX ? matrix->rowSlice[rowRepresentative[part][1] + 1] : SIZE_MAX
  };
  size_t column[3];
  for (short i = 0; i < 3; ++i)
    column[i] = entry[i] < beyond[i] ? matrix->entryColumns[entry[i]] : SIZE_MAX;

  bool compatible[4] = {
    true,
    rowRepresentative[part][0] < SIZE_MAX,
    rowRepresentative[part][1] < SIZE_MAX,
    rowRepresentative[part][1] < SIZE_MAX
  };

  for (short r = 0; r < 4; ++r)
  {
    if (compatible[r])
    {
      CMRdbgMsg(12, "Initially, %s %c%ld for part %d is compatible with %d.\n", isRow ? "Row" : "Column", isRow ? 'r' : 'c',
        row+1, part, r);
      rowData[row].type[part] = r;
    }
  }

  while (column[0] != SIZE_MAX || column[1] != SIZE_MAX || column[2] != SIZE_MAX)
  {
    size_t minColumn = (column[0] < column[1]) ? column[0] : column[1];
    if (column[2] < minColumn)
      minColumn = column[2];
    bool nonzero[3] = { column[0] == minColumn, column[1] == minColumn, column[2] == minColumn };

    if (columnData[minColumn].part == 1-part)
    {
      CMRdbgMsg(14, "column is c%ld. nonzeros are %d,%d,%d\n", minColumn+1, nonzero[0], nonzero[1], nonzero[2]);

      if (nonzero[0])
        compatible[0] = false;
      if (nonzero[0] != nonzero[1])
        compatible[1] = false;
      if (nonzero[0] != nonzero[2])
        compatible[2] = false;
      if (((nonzero[0] ? 1 : 0) + (nonzero[1] ? 1 : 0) + (nonzero[2] ? 1 : 0)) % 2)
        compatible[3] = false;
    }

    /* Advance entries with a nonzero. */
    for (short i = 0; i < 3; ++i)
    {
      if (nonzero[i])
      {
        entry[i]++;
        column[i] = entry[i] < beyond[i] ? matrix->entryColumns[entry[i]] : SIZE_MAX;
      }
    }
  }

  rowData[row].type[part] = -1;
  for (short r = 0; r < 4; ++r)
  {
    if (compatible[r])
    {
      CMRdbgMsg(12, "%s %c%ld for part %d is compatible with %d.\n", isRow ? "Row" : "Column", isRow ? 'r' : 'c',
        row+1, part, r);
      rowData[row].type[part] = r;
    }
  }

  assert((compatible[0] ? 1 : 0) + (compatible[1] ? 1 : 0) + (compatible[2] ? 1 : 0) + (compatible[3] ? 1 : 0) <= 1);
}


/**
 * \brief Assigns \p row to \p part, updating \p rowData and \p columnData.
 *
 * \returns \c true if one of the unassigned columns cannot be assigned to any part.
 */

static
bool assignRow(
  CMR_CHRMAT* matrix,                 /**< Matrix. */
  ElementData* rowData,               /**< Row element data. */
  ElementData* columnData,            /**< Column element data. */
  size_t columnRepresentative[2][2],  /**< Column representatives for each part. */
  CMR_ELEMENT* queue,                 /**< Queue. */
  size_t* pqueueBeyond,               /**< Pointer for the end of the queue. */
  size_t row,                         /**< Row to be added. */
  short part,                         /**< Part the row is to be added to. */
  bool isRow                          /**< Whether we're actually assigning a row. */
)
{
  assert(matrix);
  assert(rowData);
  assert(columnData);
  assert(queue);
  assert(row < matrix->numRows);
  assert(rowData[row].part == -1);
  assert(rowData[row].type[part] >= 0);

  bool nonzeros[4];
  nonzeros[0] = false;
  size_t entry;
  if (columnRepresentative[1-part][0] < SIZE_MAX)
  {
    CMR_CALL( CMRchrmatFindEntry(matrix, row, columnRepresentative[1-part][0], &entry) );
    nonzeros[1] = (entry < SIZE_MAX);
  }
  if (columnRepresentative[1-part][1] < SIZE_MAX)
  {
    CMR_CALL( CMRchrmatFindEntry(matrix, row, columnRepresentative[1-part][1], &entry) );
    nonzeros[2] = (entry < SIZE_MAX);
    nonzeros[3] = (nonzeros[1] != nonzeros[2]);
  }

  entry = matrix->rowSlice[row];
  size_t beyond = matrix->rowSlice[row + 1];
  size_t entryColumn = (entry < beyond) ? matrix->entryColumns[entry] : SIZE_MAX;
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    bool nonzero = column == entryColumn;

    if (columnData[column].part < 0 && columnData[column].type[1-part] >= 0)
    {
      if (nonzeros[columnData[column].type[1-part]] != nonzero)
      {
        columnData[column].type[1-part] = -1;
        if (columnData[column].type[part] >= 0)
        {
          CMRdbgMsg(14, "Adding %s %c%ld to queue.\n", isRow ? "column" : "row", isRow ? 'c' : 'r', column);
          queue[*pqueueBeyond] = isRow ? CMRcolumnToElement(column) : CMRrowToElement(column);
          (*pqueueBeyond)++;
        }
        else
        {
          CMRdbgMsg(14, "%s %c%ld cannot be assigned to any part.\n", isRow ? "Column" : "Row", isRow ? 'c' : 'r',
            column);
          return true;
        }
      }
    }

    if (nonzero)
    {
      ++entry;
      entryColumn = (entry < beyond) ? matrix->entryColumns[entry] : SIZE_MAX;
    }
  }

  rowData[row].part = part;

  return false;
}

static
CMR_ERROR extendMinorSeparation(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,       /**< Matrix. */
  CMR_CHRMAT* transpose,    /**< Transpose of \p matrix. */
  ElementData* rowData,     /**< Row element data. */
  ElementData* columnData,  /**< Column element data. */
  size_t* partRows[2],      /**< For each part, an array with its already assigned rows. */
  size_t partNumRows[2],    /**< For each part, the number of already assigned rows. */
  size_t* partColumns[2],   /**< For each part, an array with its already assigned columns. */
  size_t partNumColumns[2], /**< For each part, the number of already assigned columns. */
  CMR_ELEMENT* queue,       /**< Memory for a queue of elements. */
  CMR_SEPA** pseparation    /**< Pointer for storing the 3-separation or \c NULL if none was found. */
)
{
  assert(matrix);
  assert(transpose);
  assert(partRows);
  assert(partRows[0]);
  assert(partRows[1]);
  assert(partNumRows);
  assert(partColumns);
  assert(partColumns[0]);
  assert(partColumns[1]);
  assert(partNumColumns);
  assert(queue);
  assert(pseparation);

  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;

  for (size_t row = 0; row < numRows; ++row)
    rowData[row].part = -1;

  for (size_t column = 0; column < numColumns; ++column)
    columnData[column].part = -1;

  /* Set already asigned rows and columns.. */
  for (short part = 0; part < 2; ++part)
  {
    for (size_t r = 0; r < partNumRows[part]; ++r)
      rowData[partRows[part][r]].part = part;
    for (size_t c = 0; c < partNumColumns[part]; ++c)
      columnData[partColumns[part][c]].part = part;
  }

  size_t rowRepresentative[2][2] = { {SIZE_MAX, SIZE_MAX}, {SIZE_MAX, SIZE_MAX} };
  size_t columnRepresentative[2][2] = { {SIZE_MAX, SIZE_MAX}, {SIZE_MAX, SIZE_MAX} };

  CMRdbgMsg(12, "Checking the ranks for the following matrix:\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
  for (size_t row = 0; row < numRows; ++row)
    CMRdbgMsg(14, "Initially, row r%ld belongs to part %d.\n", row+1, rowData[row].part);
  for (size_t column = 0; column < numColumns; ++column)
    CMRdbgMsg(14, "Initially, column c%ld belongs to part %d.\n", column+1, columnData[column].part);

  size_t totalRank = 0;
  if (findRank1(matrix, rowData, columnData, rowRepresentative, columnRepresentative, 0))
  {
    CMRdbgMsg(12, "Top-right part has rank at least 1.\n");
    if (findRank2(matrix, rowData, columnData, rowRepresentative, columnRepresentative, 0))
    {
      CMRdbgMsg(12, "Top-right part has rank at least 2.\n");
      if (findRank3(matrix, rowData, columnData, rowRepresentative, 0))
      {
        CMRdbgMsg(12, "Top-right part has rank at least 3.\n");
        return CMR_OKAY;
      }
      totalRank = 2;
    }
    else
    {
      CMRdbgMsg(12, "Top-right part has rank 1.\n");
      totalRank = 1;
    }
  }
  else
  {
    CMRdbgMsg(12, "Top-right part has rank 0.\n");
  }

  if (findRank1(matrix, rowData, columnData, rowRepresentative, columnRepresentative, 1))
  {
    CMRdbgMsg(12, "Bottom-left part has rank at least 1.\n");
    ++totalRank;
    if (totalRank >= 3)
      return CMR_OKAY;
    if (findRank2(matrix, rowData, columnData, rowRepresentative, columnRepresentative, 1))
    {
      CMRdbgMsg(12, "Bottom-left part has rank at least 2.\n");
      ++totalRank;
      if (totalRank >= 3)
        return CMR_OKAY;
      if (findRank3(matrix, rowData, columnData, rowRepresentative, 1))
      {
        CMRdbgMsg(12, "Bottom-left part has rank at least 3.\n");
        return CMR_OKAY;
      }
    }
    else
    {
      CMRdbgMsg(12, "Bottom-left part has rank 1.\n");
    }
  }
  else
  {
    CMRdbgMsg(12, "Bottom-left part has rank 0.\n");
  }

  if (totalRank != 2)
  {
    CMRdbgMsg(12, "Total rank is not 2.\n");
    return CMR_OKAY;
  }

  CMRdbgMsg(12, "Total rank is 2.\n");
  CMRdbgMsg(12, "Row representatives are:    r%ld and r%ld (top-right) and r%ld and r%ld (bottom-left).\n",
    rowRepresentative[0][0]+1, rowRepresentative[0][1]+1, rowRepresentative[1][0]+1, rowRepresentative[1][1]+1);
  CMRdbgMsg(12, "Column representatives are: c%ld and c%ld (top-right) and c%ld and c%ld (bottom-left).\n",
    columnRepresentative[1][0]+1, columnRepresentative[1][1]+1, columnRepresentative[0][0]+1,
    columnRepresentative[0][1]+1);

  /* Total rank is 2. We now have to determine the types of the unassigned rows. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (rowData[row].part == -1)
    {
      determineType(matrix, rowData, columnData, rowRepresentative, row, 0, true);
      determineType(matrix, rowData, columnData, rowRepresentative, row, 1, true);
      CMRdbgMsg(12, "Row r%ld has type %d for part 0 and type %d for part 1.\n", row+1, rowData[row].type[0],
        rowData[row].type[1]);
      if (rowData[row].type[0] < 0 && rowData[row].type[1] < 0)
        return CMR_OKAY;
    }
    else
    {
      CMRdbgMsg(12, "Row r%ld belongs to part %d.\n", row+1, rowData[row].part);
    }
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    if (columnData[column].part == -1)
    {
      determineType(transpose, columnData, rowData, columnRepresentative, column, 0, false);
      determineType(transpose, columnData, rowData, columnRepresentative, column, 1, false);
      CMRdbgMsg(12, "Column c%ld has type %d for part 0 and type %d for part 1.\n", column+1,
        columnData[column].type[0], columnData[column].type[1]);
      if (columnData[column].type[0] < 0 && columnData[column].type[1] < 0)
        return CMR_OKAY;
    }
    else
    {
      CMRdbgMsg(12, "Column c%ld belongs to part %d.\n", column+1, columnData[column].part);
    }
  }

  /* Collect all elements for which there is only one part possible. */
  size_t queueFirst = 0;
  size_t queueBeyond = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (rowData[row].part >= 0)
      continue;
    if ((rowData[row].type[0] < 0 && rowData[row].type[1] >= 0)
      || (rowData[row].type[0] >= 0 && rowData[row].type[1] < 0))
    {
      queue[queueBeyond] = CMRrowToElement(row);
      ++queueBeyond;
    }
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    if (columnData[column].part >= 0)
      continue;
    if ((columnData[column].type[0] < 0 && columnData[column].type[1] >= 0)
      || (columnData[column].type[0] >= 0 && columnData[column].type[1] < 0))
    {
      queue[queueBeyond] = CMRcolumnToElement(column);
      ++queueBeyond;
    }
  }

  /* Iteratively assign elements for which there is only one part possible. */
  while (queueFirst < queueBeyond)
  {
    CMR_ELEMENT element = queue[queueFirst++];
    CMRdbgMsg(12, "Processing queue element %s.\n", CMRelementString(element, NULL));
    if (CMRelementIsRow(element))
    {
      size_t row = CMRelementToRowIndex(element);
      short part = rowData[row].type[0] >= 0 ? 0 : 1;
      assert(rowData[row].type[part] >= 0);
      assert(rowData[row].type[1-part] < 0);
      if (assignRow(matrix, rowData, columnData, columnRepresentative, queue,
        &queueBeyond, row, part, true))
      {
        return CMR_OKAY;
      }
    }
    else
    {
      size_t column = CMRelementToColumnIndex(element);
      short part = columnData[column].type[0] >= 0 ? 0 : 1;
      assert(columnData[column].type[part] >= 0);
      assert(columnData[column].type[1-part] < 0);
      if (assignRow(transpose, columnData, rowData, rowRepresentative, queue,
        &queueBeyond, column, part, false))
      {
        return CMR_OKAY;
      }
    }
  }

  /* The only unassigned elements should be compatible with both parts. So we assign them to the (potentially smaller)
   * part 0. */

  size_t countElements[2] = {0, 0};
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (rowData[row].part < 0)
    {
      CMRdbgMsg(12, "Row r%ld is unassigned and has types %d and %d. Assigning it to part 0.\n", row+1,
        rowData[row].type[0], rowData[row].type[1]);
      rowData[row].part = 0;
    }
    else
      CMRdbgMsg(12, "Row r%ld is assigned to part %d.\n", row+1, rowData[row].part);
    countElements[rowData[row].part]++;
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    if (columnData[column].part < 0)
    {
      CMRdbgMsg(12, "Column c%ld is unassigned and has types %d and %d. Assigning it to part 0.\n", column+1,
        columnData[column].type[0], columnData[column].type[1]);
      columnData[column].part = 0;
    }
    else
      CMRdbgMsg(12, "Column c%ld is assigned to part %d.\n", column+1, columnData[column].part);
    countElements[columnData[column].part]++;
  }

  CMRdbgMsg(12, "The parts of the 3-separation have %zu and %zu elements, respectively.\n", countElements[0],
    countElements[1]);
  if (countElements[0] < 4 || countElements[1] < 4)
    return CMR_OKAY;

  CMR_CALL( CMRsepaCreate(cmr, matrix->numRows, matrix->numColumns, pseparation) );
  CMR_SEPA* separation = *pseparation;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    assert(rowData[row].part == 0 || rowData[row].part == 1);
    separation->rowsFlags[row] = rowData[row].part == 0 ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    assert(columnData[column].part == 0 || columnData[column].part == 1);
    separation->columnsFlags[column] = columnData[column].part == 0 ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRregularityNestedMinorSequenceSearchThreeSeparation(CMR* cmr, DecompositionTask* task,
  DecompositionTask** punprocessed)
{
  assert(cmr);
  assert(task);
  assert(punprocessed);

  CMR_MATROID_DEC* dec = task->dec;
  assert(dec);
  assert(dec->matrix);
  assert(dec->nestedMinorsMatrix);
  assert(dec->nestedMinorsTranspose);
  assert(dec->nestedMinorsLength >= 2);

#if defined(CMR_DEBUG)
  CMRdbgMsg(6, "Searching for 3-separations along the sequence of nested minors for the following matrix:\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stdout, '0', true) );
  CMRdbgMsg(8, "Matrix with nested minor sequence:\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, dec->nestedMinorsMatrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  size_t firstNonCoGraphicMinor = (dec->nestedMinorsLastGraphic > dec->nestedMinorsLastCographic)
    ? (dec->nestedMinorsLastGraphic + 1) : (dec->nestedMinorsLastCographic + 1);
  size_t firstMinor; /* Index of first minor with at least 8 elements. */
  for (firstMinor = 1; firstMinor < dec->nestedMinorsLength; ++firstMinor)
  {
    if (dec->nestedMinorsSequenceNumRows[firstMinor] + dec->nestedMinorsSequenceNumColumns[firstMinor] >= 8)
      break;
  }

  /* We have fewer than 8 elements. */
  if (firstMinor == dec->nestedMinorsLength)
  {
    CMRdbgMsg(8, "-> irregular since fewer than 8 elements in total.\n");

    assert(false); /* TODO: add a test. */

    dec->type = CMR_MATROID_DEC_TYPE_IRREGULAR;

    /* Free the task. */
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );

    goto cleanupSequence;
  }

  /* The first minor is already neither graphic nor cographic. */
  size_t firstNonCoGraphicMinorSize = dec->nestedMinorsSequenceNumRows[firstNonCoGraphicMinor]
      + dec->nestedMinorsSequenceNumColumns[firstNonCoGraphicMinor];
  if (firstNonCoGraphicMinorSize < 9) /* TODO: this should even be 12. */
  {
    CMRdbgMsg(8, "-> irregular since minor with %zu elements is neither graphic nor cographic.\n",
      firstNonCoGraphicMinorSize);

    dec->type = CMR_MATROID_DEC_TYPE_IRREGULAR;

    /* Free the task. */
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );

    goto cleanupSequence;
  }

  CMRdbgMsg(8, "-> searching for induced 3-separations for minor indices in [%zu,%zu].\n", firstMinor,
    firstNonCoGraphicMinor);

  /* Prepare calls to the separation-extension algorithm. */
  ElementData* rowData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowData, dec->matrix->numRows) );
  ElementData* columnData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnData, dec->matrix->numColumns) );

  assert(firstMinor <= firstNonCoGraphicMinor);

  size_t* partRows[2] = { NULL, NULL };
  size_t partNumRows[2];
  size_t* partColumns[2] = { NULL, NULL };;
  size_t partNumColumns[2];
  CMR_ELEMENT* queueMemory = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &partRows[0], dec->nestedMinorsSequenceNumRows[firstNonCoGraphicMinor]) );
  CMR_CALL( CMRallocStackArray(cmr, &partRows[1], dec->nestedMinorsSequenceNumRows[firstNonCoGraphicMinor]) );
  CMR_CALL( CMRallocStackArray(cmr, &partColumns[0], dec->nestedMinorsSequenceNumColumns[firstNonCoGraphicMinor]) );
  CMR_CALL( CMRallocStackArray(cmr, &partColumns[1], dec->nestedMinorsSequenceNumColumns[firstNonCoGraphicMinor]) );
  CMR_CALL( CMRallocStackArray(cmr, &queueMemory, dec->matrix->numRows + dec->matrix->numColumns) );

  size_t firstMinorNumRows = dec->nestedMinorsSequenceNumRows[firstMinor];
  size_t firstMinorNumColumns = dec->nestedMinorsSequenceNumColumns[firstMinor];

  CMRdbgMsg(8, "Initial minor has %zu rows and %zu columns.\n", firstMinorNumRows, firstMinorNumColumns);
  CMRdbgMsg(8, "First non-(co)graphic minor has %zu rows and %zu columns.\n",
    dec->nestedMinorsSequenceNumRows[firstNonCoGraphicMinor],
    dec->nestedMinorsSequenceNumColumns[firstNonCoGraphicMinor]);

  if (task->stats)
    task->stats->enumerationCount++;

  /* Enumerate all cardinality-at-most-half subsets of the element set of the first minor. */
  short beyondBits = 1 << (firstMinorNumRows + firstMinorNumColumns);
  CMR_SEPA* separation = NULL;
  for (short bits = 0; bits < beyondBits && !separation; ++bits)
  {
    partNumRows[0] = 0;
    partNumRows[1] = 0;
    for (size_t row = 0; row < firstMinorNumRows; ++row)
    {
      short part = (bits & (1 << row)) ? 1 : 0;
      partRows[part][partNumRows[part]++] = row;
    }

    partNumColumns[0] = 0;
    partNumColumns[1] = 0;
    for (size_t column = 0; column < firstMinorNumColumns; ++column)
    {
      short part = (bits & (1 << (firstMinorNumRows + column))) ? 1 : 0;
      partColumns[part][partNumColumns[part]++] = column;
    }
    if (partNumRows[0] + partNumColumns[0] > partNumRows[1] + partNumColumns[1])
      continue;

    CMRdbgMsg(10, "Considering initial minor partition in which the first part has %zu rows and %zu columns.\n",
      partNumRows[0], partNumColumns[0]);

    if (task->stats)
      task->stats->enumerationCandidatesCount++;

    CMR_CALL( extendMinorSeparation(cmr, dec->nestedMinorsMatrix, dec->nestedMinorsTranspose, rowData, columnData,
      partRows, partNumRows, partColumns, partNumColumns, queueMemory, &separation) );
  }

  if (!separation)
  {
    /* Enumerate all subsets of elements of later minors that have at most 1 element from the previous minor and
    * at least one new. */

    for (size_t minor = firstMinor+1; minor <= firstNonCoGraphicMinor && !separation; ++minor)
    {
//       double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//       if (remainingTime < 0)
//       {
//         CMR_CALL( CMRfreeStackArray(cmr, &queueMemory) );
//         CMR_CALL( CMRfreeStackArray(cmr, &partColumns[1]) );
//         CMR_CALL( CMRfreeStackArray(cmr, &partColumns[0]) );
//         CMR_CALL( CMRfreeStackArray(cmr, &partRows[1]) );
//         CMR_CALL( CMRfreeStackArray(cmr, &partRows[0]) );
//         CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
//         CMR_CALL( CMRfreeStackArray(cmr, &rowData) );
//         return CMR_ERROR_TIMEOUT;
//       }

      CMRdbgMsg(8, "Next minor has %zu rows and %zu columns.\n", dec->nestedMinorsSequenceNumRows[minor],
        dec->nestedMinorsSequenceNumColumns[minor]);

      size_t numOldRows = dec->nestedMinorsSequenceNumRows[minor-1];
      size_t numOldColumns = dec->nestedMinorsSequenceNumColumns[minor-1];
      for (size_t old = 0; (old <= numOldRows + numOldColumns) && !separation; ++old)
      {
        partNumRows[0] = 0;
        partNumRows[1] = 0;
        partNumColumns[0] = 0;
        partNumColumns[1] = 0;

        CMR_ELEMENT oldElement = (old < numOldRows) ? CMRrowToElement(old) :
          ((old < numOldRows + numOldColumns) ? CMRcolumnToElement(old - numOldRows) : 0 );

        /* Distribute previous minors' rows to part 1 unless equal to oldElement. */
        for (size_t row = 0; row < numOldRows; ++row)
        {
          if (CMRrowToElement(row) == oldElement)
            partRows[0][partNumRows[0]++] = row;
          else
            partRows[1][partNumRows[1]++] = row;
        }

        /* Distribute previous minors' columns to part 1 unless equal to oldElement. */
        for (size_t column = 0; column < numOldColumns; ++column)
        {
          if (CMRcolumnToElement(column) == oldElement)
            partColumns[0][partNumColumns[0]++] = column;
          else
            partColumns[1][partNumColumns[1]++] = column;
        }

        CMRdbgMsg(10, "Created partition of previous minors in which the first part has %zu rows and %zu columns.\n",
          partNumRows[0], partNumColumns[0]);

        /* Enumerate distribution of new elements.
         * We store the numbers of rows and columns of the old one to skip recalculation. */

        size_t savedPartNumRows[2] = { partNumRows[0], partNumRows[1] };
        size_t savedPartNumColumns[2] = { partNumColumns[0], partNumColumns[1] };
        size_t numNewRows = dec->nestedMinorsSequenceNumRows[minor] - dec->nestedMinorsSequenceNumRows[minor-1];
        size_t numNewColumns = dec->nestedMinorsSequenceNumColumns[minor] - dec->nestedMinorsSequenceNumColumns[minor-1];
        short beyondBits = (1 << (numNewRows + numNewColumns)) - 1; /* Subtracting 1 ensures that part 0 is non-empty. */
        for (short bits = 0; bits < beyondBits && !separation; ++bits)
        {
          for (size_t newRow = 0; newRow < numNewRows; ++newRow)
          {
            size_t row = dec->nestedMinorsSequenceNumRows[minor-1] + newRow;
            short part = (bits & (1 << newRow)) ? 1 : 0;
            partRows[part][partNumRows[part]++] = row;
          }

          for (size_t newColumn = 0; newColumn < numNewColumns; ++newColumn)
          {
            size_t column = dec->nestedMinorsSequenceNumColumns[minor-1] + newColumn;
            short part = (bits & (1 << (numNewRows + newColumn))) ? 1 : 0;
            partColumns[part][partNumColumns[part]++] = column;
          }

          /* Skip this distribution if no new element chosen. */
          CMRdbgMsg(12, "Considering partition in which the first part has %zu rows and %zu columns.\n", partNumRows[0],
            partNumColumns[0]);

          if (task->stats)
            task->stats->enumerationCandidatesCount++;

          CMR_CALL( extendMinorSeparation(cmr, dec->nestedMinorsMatrix, dec->nestedMinorsTranspose, rowData, columnData,
            partRows, partNumRows, partColumns, partNumColumns, queueMemory, &separation) );

          /* Restore distribution of previous minors. */
          partNumRows[0] = savedPartNumRows[0];
          partNumRows[1] = savedPartNumRows[1];
          partNumColumns[0] = savedPartNumColumns[0];
          partNumColumns[1] = savedPartNumColumns[1];
        }
      }
    }
  }

  if (task->stats)
  {
//     task->stats->enumerationTime += (clock() - time) * 1.0 / CLOCKS_PER_SEC;
  }

  if (separation)
  {
    CMRdbgMsg(8, "Transforming 3-separation of nested minors matrix into 3-separation of node's matrix.\n");

    /* We transform the 3-separation of dec->nestedMinorsMatrix to one of dec->matrix. */

    CMR_SEPA* originalSeparation = NULL;
    CMR_CALL( CMRsepaCreate(cmr, dec->numRows, dec->numColumns, &originalSeparation) );

    for (size_t nestedRow = 0; nestedRow < dec->numRows; ++nestedRow)
    {
      CMR_ELEMENT element = dec->nestedMinorsRowsOriginal[nestedRow];
      if (CMRelementIsRow(element))
        originalSeparation->rowsFlags[CMRelementToRowIndex(element)] = separation->rowsFlags[nestedRow];
      else
        originalSeparation->columnsFlags[CMRelementToColumnIndex(element)] = separation->rowsFlags[nestedRow];
    }

    for (size_t nestedColumn = 0; nestedColumn < dec->numColumns; ++nestedColumn)
    {
      CMR_ELEMENT element = dec->nestedMinorsColumnsOriginal[nestedColumn];
      if (CMRelementIsRow(element))
        originalSeparation->rowsFlags[CMRelementToRowIndex(element)] = separation->columnsFlags[nestedColumn];
      else
        originalSeparation->columnsFlags[CMRelementToColumnIndex(element)] = separation->columnsFlags[nestedColumn];
    }

    CMR_CALL( CMRsepaFree(cmr, &separation) );

    if (dec->transpose == NULL)
    {
      CMR_CALL( CMRchrmatTranspose(cmr, dec->matrix, &dec->transpose) );
    }
    CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, originalSeparation, dec->matrix, dec->transpose, NULL) );

    if (dec->isTernary)
    {
      bool isTernary;
      CMR_SUBMAT* violatorSubmatrix = NULL;
      CMR_CALL( CMRsepaCheckTernary(cmr, originalSeparation, dec->matrix, &isTernary, &violatorSubmatrix) );

      if (!isTernary)
      {
        CMRdbgMsg(8, "-> 2x2 submatrix with bad determinant.\n");

        CMR_CALL( CMRmatroiddecUpdateSubmatrix(cmr, dec, violatorSubmatrix, CMR_MATROID_DEC_TYPE_DETERMINANT) );
        assert(dec->type != CMR_MATROID_DEC_TYPE_DETERMINANT);

        CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
        CMR_CALL( CMRsepaFree(cmr, &originalSeparation) );

        goto cleanupSearch;
      }
    }

    CMR_CALL( CMRregularityDecomposeThreeSum(cmr, task, punprocessed, originalSeparation) );

    CMR_CALL( CMRsepaFree(cmr, &originalSeparation) );
  }
  else
  {
    assert(!"NOT TESTED"); // TODO: Add a unittest.
    task->dec->type = CMR_MATROID_DEC_TYPE_IRREGULAR;

    /* Free the task. */
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }

cleanupSearch:

  /* Free working storage. */
  CMR_CALL( CMRfreeStackArray(cmr, &queueMemory) );
  CMR_CALL( CMRfreeStackArray(cmr, &partColumns[1]) );
  CMR_CALL( CMRfreeStackArray(cmr, &partColumns[0]) );
  CMR_CALL( CMRfreeStackArray(cmr, &partRows[1]) );
  CMR_CALL( CMRfreeStackArray(cmr, &partRows[0]) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowData) );

cleanupSequence:

  /* Free the sequence of nested minors. */
  CMR_CALL( CMRchrmatFree(cmr, &dec->nestedMinorsMatrix) );
  CMR_CALL( CMRchrmatFree(cmr, &dec->nestedMinorsTranspose) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsSequenceNumRows) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsSequenceNumColumns) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsRowsOriginal) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsColumnsOriginal) );

  return CMR_OKAY;
}

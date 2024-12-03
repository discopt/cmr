// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/separation.h>

#include "env_internal.h"
#include "bipartite_graph.h"

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

CMR_ERROR CMRsepaCreate(CMR* cmr, size_t numRows, size_t numColumns, CMR_SEPA** psepa)
{
  assert(cmr);
  assert(psepa);

  CMR_CALL( CMRallocBlock(cmr, psepa) );
  CMR_SEPA* sepa = *psepa;
  sepa->type = 0;
  sepa->numRows = numRows;
  sepa->numColumns = numColumns;
  sepa->rowsFlags = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &sepa->rowsFlags, numRows) );
  sepa->columnsFlags = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &sepa->columnsFlags, numColumns) );

  return CMR_OKAY;
}

CMR_ERROR CMRsepaFree(CMR* cmr, CMR_SEPA** psepa)
{
  assert(cmr);
  assert(psepa);

  if (*psepa)
  {
    CMR_SEPA* sepa = *psepa;
    CMR_CALL( CMRfreeBlockArray(cmr, &sepa->rowsFlags) );
    CMR_CALL( CMRfreeBlockArray(cmr, &sepa->columnsFlags) );
    CMR_CALL( CMRfreeBlock(cmr, psepa) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRsepaComputeSizes(CMR_SEPA* sepa, size_t* pnumRowsTopLeft, size_t* pnumColumnsTopLeft,
  size_t* pnumRowsBottomRight, size_t* pnumColumnsBottomRight)
{
  assert(sepa);

  size_t numRowsTopLeft = 0;
  size_t numColumnsTopLeft = 0;
  size_t numRowsBottomRight = 0;
  size_t numColumnsBottomRight = 0;

  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    if ((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      ++numRowsTopLeft;
    else
      ++numRowsBottomRight;
  }

  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    if ((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      ++numColumnsTopLeft;
    else
      ++numColumnsBottomRight;
  }

  if (pnumRowsTopLeft)
    *pnumRowsTopLeft = numRowsTopLeft;
  if (pnumColumnsTopLeft)
    *pnumColumnsTopLeft = numColumnsTopLeft;
  if (pnumRowsBottomRight)
    *pnumRowsBottomRight = numRowsBottomRight;
  if (pnumColumnsBottomRight)
    *pnumColumnsBottomRight = numColumnsBottomRight;

  return CMR_OKAY;
}

typedef struct
{
  uint8_t considered;
  int8_t reprNonzero[8]; /* r_1, -r_1, r_2, -r_2, r_1+r_2, -r_1-r_2, r_1-r_2, -r_1+r_2 */
} ElementData;

/**
 * \brief Compute the binary and ternary of a submatrix of \p matrix (if at most 2).
 */

static
CMR_ERROR computeSubmatrixRank(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Matrix. */
  size_t numRows,                 /**< Number of rows. */
  size_t* rows,                   /**< Array with row indices. */
  bool* columnsConsidered,        /**< Array indicating for each column of \p matrix whether it shall be considered. */
  CMR_SEPA_FLAGS* rowsFlags,      /**< Array of length \p numRows for storing the flags for representatives. */
  size_t* prank,                  /**< Pointer for storing the computed rank. */
  CMR_SUBMAT** pviolatorSubmatrix /**< Pointer for storing a submatrix of absolute determinant 2 (or \c NULL). */
)
{
  assert(matrix);
  assert(rows);
  assert(columnsConsidered);
  assert(rowsFlags);

  size_t rank = 0;
  size_t theRows[3]; /* The (up to two) representative rows and the current row. */
  size_t theColumns[3] = { SIZE_MAX, SIZE_MAX, SIZE_MAX }; /* The (up to two) representative columns and the current column. */
  size_t beyond[3];
  size_t entry[3];
  size_t column[3];
  int8_t values[3] = { INT8_MAX, INT8_MAX, INT8_MAX };

  CMRdbgMsg(12, "computeSubmatrixRank() for a submatrix with %zu rows.\n", numRows);

  for (size_t r = 0; r < numRows && rank <= 2; ++r)
  {
    theRows[rank] = rows[r];

    CMRdbgMsg(14, "Considering row r%zu.\n", rows[r]+1);

    /* Prepare data for iterating simultaneously through all relevant rows. */
    for (size_t i = 0; i <= rank; ++i)
    {
      beyond[i] = matrix->rowSlice[theRows[i] + 1];
      entry[i] = matrix->rowSlice[theRows[i]];
      column[i] = entry[i] < beyond[i] ? matrix->entryColumns[entry[i]] : SIZE_MAX;
    }

    /* Indicators (indexed by the representative vectors). */
    bool zero = true;
    bool binary_x = rank >= 1;
    bool binary_y = rank >= 2;
    bool binary_x_plus_y = rank >= 2;

    bool ternary_x = rank >= 1;
    bool ternary_minus_x = rank >= 1;
    bool ternary_y = rank >= 2;
    size_t ternary_y_column = SIZE_MAX;
    bool ternary_minus_y = rank >= 2;
    bool ternary_x_plus_y = rank >= 2;
    bool ternary_x_minus_y = rank >= 2;
    bool ternary_minus_x_plus_y = rank >= 2;
    bool ternary_minus_x_minus_y = rank >= 2;

    CMRdbgMsg(16, "Rank = %zu; initial indicators are %d|%d|%d (binary) and %d%d|%d%d|%d%d%d%d (ternary).\n", rank,
      binary_x, binary_y, binary_x_plus_y, ternary_x, ternary_minus_x, ternary_y, ternary_minus_y,
      ternary_x_plus_y, ternary_x_minus_y, ternary_minus_x_plus_y, ternary_minus_x_minus_y);

    /* Loop over all columns of the representative rows. */
    while ( (rank > 0 && column[0] < SIZE_MAX) || (rank > 1 && column[1] < SIZE_MAX) || (column[rank] < SIZE_MAX) )
    {
      CMRdbgMsg(16, "Iteration:\n");
      for (size_t i = 0; i <= rank; ++i)
      {
        CMRdbgMsg(18, "Entry #%zu is #%zu in [#%zu,#%zu) (c%zu).\n", i, entry[i], matrix->rowSlice[theRows[i]],
          beyond[i], column[i]+1);
      }

      /* Which column are we looking at? */
      size_t currentColumn = SIZE_MAX;
      for (size_t i = 0; i <= rank; ++i)
      {
        if (entry[i] < beyond[i] && column[i] < currentColumn)
          currentColumn = column[i];
      }
      if (currentColumn == SIZE_MAX)
        break;

      if (columnsConsidered[currentColumn])
      {
        for (size_t i = 0; i <= rank; ++i)
          values[i] = (column[i] == currentColumn) ? matrix->entryValues[entry[i]] : 0;

        if (values[rank])
          zero = false;

        /* Update indicators for each vector in the binary span of the found ones. */
        binary_x = binary_x && (abs(values[0]) - abs(values[rank]) == 0);
        binary_y = binary_y && (abs(values[1]) - abs(values[rank]) == 0);
        binary_x_plus_y = binary_x_plus_y && ((abs(values[0]) + abs(values[1]) - abs(values[rank])) % 2 == 0);

        /* Update indicators for each vector in the ternary span of the found ones. */
        ternary_x = ternary_x && ((values[0] - values[rank]) % 3 == 0);
        ternary_minus_x = ternary_minus_x && ((-values[0] - values[rank]) % 3 == 0);
        ternary_y = ternary_y && ((values[1] - values[rank]) % 3 == 0);
        ternary_minus_y = ternary_minus_y && ((-values[1] - values[rank]) % 3 == 0);
        ternary_x_plus_y = ternary_x_plus_y && ((values[0] + values[1] - values[rank]) % 3 == 0);
        ternary_x_minus_y = ternary_x_minus_y && ((values[0] - values[1] - values[rank]) % 3 == 0);
        ternary_minus_x_plus_y = ternary_minus_x_plus_y && ((-values[0] + values[1] - values[rank]) % 3 == 0);
        ternary_minus_x_minus_y = ternary_minus_x_minus_y && ((-values[0] - values[1] - values[rank]) % 3 == 0);

        CMRdbgMsg(16, "Considering column c%zu (values %d,%d,%d); indicators are %d|%d|%d (binary) "
          "and %d%d|%d%d|%d%d%d%d (ternary).\n", currentColumn+1, values[0], values[1], values[2],
          binary_x, binary_y, binary_x_plus_y, ternary_x, ternary_minus_x, ternary_y, ternary_minus_y,
          ternary_x_plus_y, ternary_x_minus_y, ternary_minus_x_plus_y, ternary_minus_x_minus_y);

        /** Update indicator columns. */
        if ( (ternary_y_column == SIZE_MAX) && binary_y && (!ternary_y || !ternary_minus_y) )
          ternary_y_column = currentColumn;

        if (pviolatorSubmatrix)
        {
          if (binary_x && !ternary_x && !ternary_minus_x)
          {
            CMR_CALL( CMRsubmatCreate(cmr, 2, 2, pviolatorSubmatrix) );
            CMR_SUBMAT* violatorSubmatrix = *pviolatorSubmatrix;
            violatorSubmatrix->rows[0] = theRows[0];
            violatorSubmatrix->rows[1] = theRows[rank];
            violatorSubmatrix->columns[0] = theColumns[0];
            violatorSubmatrix->columns[1] = currentColumn;

            CMRdbgMsg(16, "Found 2-by-2 submatrix with absolute determinant 2: r%zu,r%zu,c%zu,c%zu.\n",
              theRows[0]+1, theRows[rank]+1, theColumns[0]+1, currentColumn+1);

            *prank = rank;
            return CMR_OKAY;
          }

          if (binary_y && !ternary_y && !ternary_minus_y)
          {
            CMR_CALL( CMRsubmatCreate(cmr, 2, 2, pviolatorSubmatrix) );
            CMR_SUBMAT* violatorSubmatrix = *pviolatorSubmatrix;
            violatorSubmatrix->rows[0] = theRows[1];
            violatorSubmatrix->rows[1] = theRows[rank];
            violatorSubmatrix->columns[0] = ternary_y_column;
            violatorSubmatrix->columns[1] = currentColumn;

            CMRdbgMsg(16, "Found 2-by-2 submatrix with absolute determinant 2: r%zu,r%zu,c%zu,c%zu.\n",
              theRows[1]+1, theRows[rank]+1, ternary_y_column+1, currentColumn+1);

            assert(ternary_y_column != currentColumn);

            *prank = rank;
            return CMR_OKAY;
          }

          if (binary_x_plus_y && !ternary_x_plus_y && !ternary_x_minus_y && !ternary_minus_x_plus_y
            && !ternary_minus_x_minus_y)
          {
            CMR_CALL( CMRsubmatCreate(cmr, 3, 3, pviolatorSubmatrix) );
            CMR_SUBMAT* violatorSubmatrix = *pviolatorSubmatrix;
            violatorSubmatrix->rows[0] = theRows[0];
            violatorSubmatrix->rows[1] = theRows[1];
            violatorSubmatrix->rows[2] = theRows[2];
            violatorSubmatrix->columns[0] = theColumns[0];
            violatorSubmatrix->columns[1] = theColumns[1];
            violatorSubmatrix->columns[2] = currentColumn;

            CMRdbgMsg(16, "Found 3-by-3 submatrix with absolute determinant 2: r%zu,r%zu,r%zu,c%zu,c%zu,c%zu.\n",
              theRows[0]+1, theRows[1]+1, theRows[2]+1, theColumns[0]+1, theColumns[1]+1, currentColumn+1);

            return CMR_OKAY;
          }
        }

        /* No need to continue if it increases the rank. */
        if (!zero && !binary_x && !binary_y && !binary_x_plus_y)
        {
          theColumns[rank] = currentColumn;
          break;
        }
      }

      /* Advance all entries whose column is the current one. */
      CMRdbgMsg(16, "Advancing subset of the entries.\n");
      for (size_t i = 0; i <= rank; ++i)
      {
        if (entry[i] < beyond[i] && column[i] == currentColumn)
        {
          entry[i]++;
          column[i] = (entry[i] < beyond[i]) ? matrix->entryColumns[entry[i]] : SIZE_MAX;
        }
      }
    }

    CMRdbgMsg(16, "Completed iteration over rows.\n");

    if (zero)
    {
      CMRdbgMsg(14, "-> zero vector.\n");
    }
    else if (binary_x)
    {
      rowsFlags[r] = CMR_SEPA_FLAG_RANK1;
      CMRdbgMsg(14, "-> 1st representative vector.\n");
    }
    else if (binary_y)
    {
      rowsFlags[r] = CMR_SEPA_FLAG_RANK2;
      CMRdbgMsg(14, "-> 2nd representative vector.\n");
    }
    else if (binary_x_plus_y)
    {
      rowsFlags[r] = CMR_SEPA_FLAG_RANK1 | CMR_SEPA_FLAG_RANK2;
      CMRdbgMsg(14, "-> 1st + 2nd representative vector.\n");
    }
    else
    {
      ++rank;
      rowsFlags[r] = (rank == 1) ? CMR_SEPA_FLAG_RANK1 : CMR_SEPA_FLAG_RANK2;
      CMRdbgMsg(14, "-> New %s representative vector.\n", rank == 1 ? "1st" : "2nd");
    }
  }

  *prank = rank;
  return CMR_OKAY;
}

/**
 * \brief Implementation of \ref CMRsepaFindBinaryRepresentatives and \ref CMRsepaFindBinaryRepresentativesSubmatrix.
 */

static
CMR_ERROR findRepresentatives(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_SEPA* sepa,                 /**< Separation. */
  CMR_CHRMAT* matrix,             /**< Matrix. */
  CMR_CHRMAT* transpose,          /**< Transpose of \p matrix. */
  size_t* submatrixRows,          /**< Array mapping a submatrix row to a row of \p matrix. */
  size_t* submatrixColumns,       /**< Array mapping a submatrix column to a column of \p matrix. */
  bool* pswapped,                 /**< Pointer for storing whether parts were swapped (may be \c NULL). */
  CMR_SUBMAT** pviolatorSubmatrix /**< Pointer for storing a violator submatrix if the ternary rank differs
                                   **< (may be \c NULL). */
)
{
  assert(cmr);
  assert(sepa);
  assert(matrix);
  assert(transpose);
  assert(submatrixRows);
  assert(sepa->numRows <= matrix->numRows);
  assert(sepa->numColumns <= matrix->numColumns);

#if defined(CMR_DEBUG)
  CMRdbgMsg(8, "Finding representatives of low-rank submatrices.\n");
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
  for (size_t r = 0; r < sepa->numRows; ++r)
  {
    size_t row = submatrixRows[r];
    CMRdbgMsg(10, "Submatrix row %zu is row r%zu belongs to part %d.\n", r, row+1,
      ((sepa->rowsFlags[r] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST) ? 0 : 1 );
  }
  for (size_t c = 0; c < sepa->numColumns; ++c)
  {
    size_t column = submatrixColumns[c];
    CMRdbgMsg(10, "Submatrix column %zu is column c%zu belongs to part %d.\n", c, column+1,
      ((sepa->columnsFlags[c] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST) ? 0 : 1 );
  }
#endif /* CMR_DEBUG */

  /* We first create the inverse mappings of submatrixRows and submatrixColumns. */
  size_t* rowsToSubmatrixRow = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowsToSubmatrixRow, matrix->numRows) );
  if (sepa->numRows < matrix->numRows)
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
      rowsToSubmatrixRow[row] = SIZE_MAX;
  }
  for (size_t r = 0; r < sepa->numRows; ++r)
  {
    rowsToSubmatrixRow[submatrixRows[r]] = r;
    sepa->rowsFlags[r] &= CMR_SEPA_MASK_CHILD; /* Clear the representative-flags. */
  }

  size_t* columnsToSubmatrixColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnsToSubmatrixColumn, matrix->numColumns) );
  if (sepa->numColumns < matrix->numColumns)
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columnsToSubmatrixColumn[column] = SIZE_MAX;
  }
  for (size_t c = 0; c < sepa->numColumns; ++c)
  {
    columnsToSubmatrixColumn[submatrixColumns[c]] = c;
    sepa->columnsFlags[c] &= CMR_SEPA_MASK_CHILD; /* Clear the representative-flags. */
  }

  /* Initialize actual algorithm. */

  size_t numMax = matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns;
  size_t numMajors;
  size_t* majors = NULL;
  uint8_t* majorsRepresentatives = NULL;
  ElementData* minorsData = NULL;
  CMR_SEPA_FLAGS* majorsFlags = NULL;
  bool* minorsConsidered = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &majors, numMax) );
  CMR_CALL( CMRallocStackArray(cmr, &majorsRepresentatives, numMax) );
  CMR_CALL( CMRallocStackArray(cmr, &minorsData, numMax) );
  CMR_CALL( CMRallocStackArray(cmr, &majorsFlags, numMax) );
  CMR_CALL( CMRallocStackArray(cmr, &minorsConsidered, numMax) );

  /* Compute bottom-left row rank. */
  numMajors = 0;
  for (size_t submatrixRow = 0; submatrixRow < sepa->numRows; ++submatrixRow)
  {
    size_t row = submatrixRows[submatrixRow];
    if ((sepa->rowsFlags[submatrixRow] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
    {
      CMRdbgMsg(12, "Major #%zu is r%zu.\n", numMajors, row+1);
      majors[numMajors] = row;
      majorsFlags[numMajors] = 0;
      ++numMajors;
    }
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    size_t submatrixColumn = columnsToSubmatrixColumn[column];
    minorsConsidered[column] = (submatrixColumn < SIZE_MAX)
      && ((sepa->columnsFlags[submatrixColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST);
  }

  size_t rankBottomLeft;
  CMRdbgMsg(10, "Computing bottom-left row-rank...\n");
  CMR_CALL( computeSubmatrixRank(cmr, matrix, numMajors, majors, minorsConsidered, majorsFlags, &rankBottomLeft,
    pviolatorSubmatrix) );

  if (pviolatorSubmatrix && *pviolatorSubmatrix)
  {
    CMR_SUBMAT* violatorSubmatrix = *pviolatorSubmatrix;

    /* Return submatrix as one of the given submatrix. */
    for (size_t row = 0; row < violatorSubmatrix->numRows; ++row)
    {
      assert( rowsToSubmatrixRow[violatorSubmatrix->rows[row]] != SIZE_MAX );
      violatorSubmatrix->rows[row] = rowsToSubmatrixRow[violatorSubmatrix->rows[row]];
    }
    for (size_t column = 0; column < violatorSubmatrix->numColumns; ++column)
    {
      assert( columnsToSubmatrixColumn[violatorSubmatrix->columns[column]] != SIZE_MAX );
      violatorSubmatrix->columns[column] = columnsToSubmatrixColumn[violatorSubmatrix->columns[column]];
    }

    goto cleanup;
  }

  CMRdbgMsg(10, "Bottom-left rank is %zu.\n", rankBottomLeft);

  /* Copy back flags to sepa. */
  numMajors = 0;
  for (size_t submatrixRow = 0; submatrixRow < sepa->numRows; ++submatrixRow)
  {
    if ((sepa->rowsFlags[submatrixRow] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
    {
      CMRdbgMsg(12, "Copy back: submatrix row %zu is major %zu with flags %d.\n", submatrixRow, numMajors,
        majorsFlags[numMajors]);
      sepa->rowsFlags[submatrixRow] |= majorsFlags[numMajors++];
    }
  }


  /* Compute top-right rank. */
  numMajors = 0;
  for (size_t submatrixRow = 0; submatrixRow < sepa->numRows; ++submatrixRow)
  {
    size_t row = submatrixRows[submatrixRow];
    if ((sepa->rowsFlags[submatrixRow] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
    {
      majors[numMajors] = row;
      majorsFlags[numMajors] = 0;
      ++numMajors;
    }
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    size_t submatrixColumn = columnsToSubmatrixColumn[column];
    minorsConsidered[column] = (submatrixColumn < SIZE_MAX)
      && ((sepa->columnsFlags[submatrixColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND);
  }

  size_t rankTopRight;
  CMRdbgMsg(10, "Computing top-right row-rank...\n");
  CMR_CALL( computeSubmatrixRank(cmr, matrix, numMajors, majors, minorsConsidered, majorsFlags, &rankTopRight,
    pviolatorSubmatrix) );

  if (pviolatorSubmatrix && *pviolatorSubmatrix)
  {
    CMR_SUBMAT* violatorSubmatrix = *pviolatorSubmatrix;

    /* Return submatrix as one of the given submatrix. */
    for (size_t row = 0; row < violatorSubmatrix->numRows; ++row)
    {
      assert( rowsToSubmatrixRow[violatorSubmatrix->rows[row]] != SIZE_MAX );
      violatorSubmatrix->rows[row] = rowsToSubmatrixRow[violatorSubmatrix->rows[row]];
    }
    for (size_t column = 0; column < violatorSubmatrix->numColumns; ++column)
    {
      assert( columnsToSubmatrixColumn[violatorSubmatrix->columns[column]] != SIZE_MAX );
      violatorSubmatrix->columns[column] = columnsToSubmatrixColumn[violatorSubmatrix->columns[column]];
    }

    goto cleanup;
  }

  CMRdbgMsg(10, "Top-right binary rank is %zu.\n", rankTopRight);

  /* Copy back flags to sepa. */
  numMajors = 0;
  for (size_t submatrixRow = 0; submatrixRow < sepa->numRows; ++submatrixRow)
  {
    if ((sepa->rowsFlags[submatrixRow] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      sepa->rowsFlags[submatrixRow] |= majorsFlags[numMajors++];
  }


  /* Compute bottom-left column rank. */
  numMajors = 0;
  for (size_t submatrixColumn = 0; submatrixColumn < sepa->numColumns; ++submatrixColumn)
  {
    size_t column = submatrixColumns[submatrixColumn];
    if ((sepa->columnsFlags[submatrixColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
    {
      CMRdbgMsg(12, "Major #%zu is c%zu.\n", numMajors, column+1);
      majors[numMajors] = column;
      majorsFlags[numMajors] = 0;
      ++numMajors;
    }
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t submatrixRow = rowsToSubmatrixRow[row];
    minorsConsidered[row] = (submatrixRow < SIZE_MAX)
      && ((sepa->rowsFlags[submatrixRow] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND);
  }

  size_t rank;
  CMRdbgMsg(10, "Computing bottom-left column-rank...\n");
  CMR_CALL( computeSubmatrixRank(cmr, transpose, numMajors, majors, minorsConsidered, majorsFlags, &rank, NULL) );
  CMRdbgMsg(10, "Bottom-left rank is confirmed to be %zu.\n", rank);
  assert(rank == rankBottomLeft);

  /* Copy back flags to sepa. */
  numMajors = 0;
  for (size_t submatrixColumn = 0; submatrixColumn < sepa->numColumns; ++submatrixColumn)
  {
    if ((sepa->columnsFlags[submatrixColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      sepa->columnsFlags[submatrixColumn] |= majorsFlags[numMajors++];
  }

  /* Compute top-right column rank. */
  numMajors = 0;
  for (size_t submatrixColumn = 0; submatrixColumn < sepa->numColumns; ++submatrixColumn)
  {
    size_t column = submatrixColumns[submatrixColumn];
    if ((sepa->columnsFlags[submatrixColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
    {
      CMRdbgMsg(12, "Major #%zu is c%zu.\n", numMajors, column+1);
      majors[numMajors] = column;
      majorsFlags[numMajors] = 0;
      ++numMajors;
    }
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t submatrixRow = rowsToSubmatrixRow[row];
    minorsConsidered[row] = (submatrixRow < SIZE_MAX)
      && ((sepa->rowsFlags[submatrixRow] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST);
  }

  CMRdbgMsg(10, "Computing top-right column-rank...\n");
  CMR_CALL( computeSubmatrixRank(cmr, transpose, numMajors, majors, minorsConsidered, majorsFlags, &rank, NULL) );
  CMRdbgMsg(10, "Top-right rank is confirmed to be %zu.\n", rank);
  assert(rank == rankTopRight);

  /* Copy back flags to sepa. */
  numMajors = 0;
  for (size_t submatrixColumn = 0; submatrixColumn < sepa->numColumns; ++submatrixColumn)
  {
    if ((sepa->columnsFlags[submatrixColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
      sepa->columnsFlags[submatrixColumn] |= majorsFlags[numMajors++];
  }

  /* Potentially swap meanings. */

  if (rankBottomLeft < rankTopRight)
  {
    /* Swap parts. */
    for (size_t row = 0; row < sepa->numRows; ++row)
    {
      sepa->rowsFlags[row] = (sepa->rowsFlags[row] & CMR_SEPA_MASK_EXTRA) |
        (((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST) ? CMR_SEPA_SECOND : CMR_SEPA_FIRST);
    }
    for (size_t column = 0; column < sepa->numColumns; ++column)
    {
      sepa->columnsFlags[column] = (sepa->columnsFlags[column] & CMR_SEPA_MASK_EXTRA) |
        (((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST) ? CMR_SEPA_SECOND : CMR_SEPA_FIRST);
    }

    if (pswapped)
      *pswapped = true;
  }
  else if (pswapped)
    *pswapped = false;

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &minorsConsidered) );
  CMR_CALL( CMRfreeStackArray(cmr, &majorsFlags) );
  CMR_CALL( CMRfreeStackArray(cmr, &minorsData) );
  CMR_CALL( CMRfreeStackArray(cmr, &majorsRepresentatives) );
  CMR_CALL( CMRfreeStackArray(cmr, &majors) );

  CMR_CALL( CMRfreeStackArray(cmr, &columnsToSubmatrixColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowsToSubmatrixRow) );

  /* Set the type. */

  if (!pviolatorSubmatrix || !*pviolatorSubmatrix)
  {
    int rank = rankBottomLeft + rankTopRight;
    if (rank == 1)
      sepa->type = CMR_SEPA_TYPE_TWO;
    else if (rankBottomLeft == 1 && rankTopRight == 1)
      sepa->type = CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS;
    else if (rank == 2)
      sepa->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRsepaFindBinaryRepresentatives(CMR* cmr, CMR_SEPA* sepa, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose,
  bool* pswapped, CMR_SUBMAT** pviolator)
{
  assert(cmr);
  assert(sepa);
  assert(matrix);

  size_t* submatrixRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &submatrixRows, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    submatrixRows[row] = row;

  size_t* submatrixColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &submatrixColumns, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    submatrixColumns[column] = column;

  CMR_CALL( findRepresentatives(cmr, sepa, matrix, transpose, submatrixRows, submatrixColumns, pswapped, pviolator) );

  CMR_CALL( CMRfreeStackArray(cmr, &submatrixColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &submatrixRows) );

  return CMR_OKAY;
}

CMR_ERROR CMRsepaFindBinaryRepresentativesSubmatrix(CMR* cmr, CMR_SEPA* sepa, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose,
  CMR_SUBMAT* submatrix, bool* pswapped, CMR_SUBMAT** pviolator)
{
  assert(cmr);
  assert(sepa);
  assert(matrix);
  assert(submatrix);

  assert(sepa->numRows <= matrix->numRows);
  assert(sepa->numColumns <= matrix->numColumns);
  assert(sepa->numRows == submatrix->numRows);
  assert(sepa->numColumns == submatrix->numColumns);

  CMR_CALL( findRepresentatives(cmr, sepa, matrix, transpose, submatrix->rows, submatrix->columns, pswapped,
    pviolator) );

  return CMR_OKAY;
}

CMR_ERROR CMRsepaGetRepresentatives(CMR_SEPA* sepa, size_t reprRows[2][3], size_t reprColumns[2][3])
{
  assert(sepa);
  assert(reprRows);
  assert(reprColumns);

  for (size_t i = 0; i < 6; ++i)
  {
    reprRows[i / 3][i % 3] = SIZE_MAX;
    reprColumns[i / 3][i % 3] = SIZE_MAX;
  }

  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    CMR_SEPA_FLAGS flags = sepa->rowsFlags[row];
    if ((flags & CMR_SEPA_MASK_EXTRA) == 0)
      continue;

    size_t child = (flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST ? 1 : 0;
    flags &= CMR_SEPA_MASK_EXTRA;
    if (flags == CMR_SEPA_FLAG_RANK1)
      reprRows[child][0] = row;
    else if (flags == CMR_SEPA_FLAG_RANK2)
      reprRows[child][1] = row;
    else if (flags == (CMR_SEPA_FLAG_RANK1 | CMR_SEPA_FLAG_RANK2))
      reprRows[child][2] = row;
  }

  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    CMR_SEPA_FLAGS flags = sepa->columnsFlags[column];
    if ((flags & CMR_SEPA_MASK_EXTRA) == 0)
      continue;

    size_t child = (flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST ? 1 : 0;
    flags &= CMR_SEPA_MASK_EXTRA;
    if (flags == CMR_SEPA_FLAG_RANK1)
      reprColumns[child][0] = column;
    else if (flags == CMR_SEPA_FLAG_RANK2)
      reprColumns[child][1] = column;
    else if (flags == (CMR_SEPA_FLAG_RANK1 | CMR_SEPA_FLAG_RANK2))
      reprColumns[child][2] = column;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRsepaGetProjection(CMR_SEPA* sepa, size_t part, size_t* rowsToPart, size_t* columnsToPart,
  size_t* pnumPartRows, size_t* pnumPartColumns)
{
  assert(sepa);
  assert(part <= 1);
  assert(rowsToPart);
  assert(columnsToPart);
  assert(pnumPartRows);
  assert(pnumPartColumns);

  CMR_SEPA_FLAGS thisPart = part ? CMR_SEPA_SECOND : CMR_SEPA_FIRST;

  size_t numPartRows = 0;
  for (size_t row = 0; row < sepa->numRows; ++row)
    rowsToPart[row] = SIZE_MAX;
  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    CMR_SEPA_FLAGS flags = sepa->rowsFlags[row];
    if ((flags & CMR_SEPA_MASK_CHILD) == thisPart)
      rowsToPart[row] = numPartRows++;
  }
  size_t previous[3] = { SIZE_MAX, SIZE_MAX, SIZE_MAX };
  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    CMR_SEPA_FLAGS flags = sepa->rowsFlags[row];
    if ((flags & CMR_SEPA_MASK_CHILD) != thisPart)
    {
      flags &= CMR_SEPA_MASK_EXTRA;
      size_t repr = (flags & CMR_SEPA_FLAG_RANK1)
        ? ((flags & CMR_SEPA_FLAG_RANK2) ? 2 : 0 )
        : ((flags & CMR_SEPA_FLAG_RANK2) ? 1 : SIZE_MAX );
      if (repr == SIZE_MAX)
        continue;

      rowsToPart[row] = numPartRows + repr;
      if (previous[repr] != SIZE_MAX)
        rowsToPart[previous[repr]] = SIZE_MAX;
      previous[repr] = row;
    }
  }
  *pnumPartRows = numPartRows;

  size_t numPartColumns = 0;
  for (size_t column = 0; column < sepa->numColumns; ++column)
    columnsToPart[column] = SIZE_MAX;
  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    CMR_SEPA_FLAGS flags = sepa->columnsFlags[column];
    if ((flags & CMR_SEPA_MASK_CHILD) == thisPart)
      columnsToPart[column] = numPartColumns++;
  }
  previous[0] = SIZE_MAX;
  previous[1] = SIZE_MAX;
  previous[2] = SIZE_MAX;
  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    CMR_SEPA_FLAGS flags = sepa->columnsFlags[column];
    if ((flags & CMR_SEPA_MASK_CHILD) != thisPart)
    {
      flags &= CMR_SEPA_MASK_EXTRA;
      size_t repr = (flags & CMR_SEPA_FLAG_RANK1)
        ? ((flags & CMR_SEPA_FLAG_RANK2) ? 2 : 0 )
        : ((flags & CMR_SEPA_FLAG_RANK2) ? 1 : SIZE_MAX );
      if (repr == SIZE_MAX)
        continue;

      columnsToPart[column] = numPartColumns + repr;
      if (previous[repr] != SIZE_MAX)
        columnsToPart[previous[repr]] = SIZE_MAX;
      previous[repr] = column;
    }
  }
  *pnumPartColumns = numPartColumns;

  return CMR_OKAY;
}

/**
 * \brief Implementation of \p CMRsepaCheckTernary and \p CMRsepaCheckTernarySubmatrix.
 */

static
CMR_ERROR checkTernary(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_SEPA* sepa,                   /**< Separation. */
  CMR_CHRMAT* matrix,               /**< Matrix. */
  size_t* submatrixRows,            /**< Array mapping a submatrix row to a row of \p matrix. */
  size_t* columnsToSubmatrixColumn, /**< Array mapping a column of \p matrix to a column of the submatrix or \c NULL. */
  bool* pisTernary,                 /**< Pointer for storing whether the check passed. */
  CMR_SUBMAT** pviolator            /**< Pointer for storing a violator submatrix (may be \c NULL). */
)
{
  assert(cmr);
  assert(sepa);
  assert(matrix);
  assert(pisTernary);

  size_t reprRows[2][3];
  size_t reprColumns[2][3];
  CMR_CALL( CMRsepaGetRepresentatives(sepa, reprRows, reprColumns) );

  if (sepa->type == CMR_SEPA_TYPE_TWO)
  {
    CMRdbgMsg(4, "Checking 2-separation for being ternary.\n");

    *pisTernary = true;

    /* We find a row belonging to the rank-1 submatrix. */
    size_t representativeSubmatrixRow = SIZE_MAX;
    for (size_t submatrixRow = 0; submatrixRow < sepa->numRows; ++submatrixRow)
    {
      CMR_SEPA_FLAGS flags = sepa->rowsFlags[submatrixRow];
      /* Skip top-left rows. */
      if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
        continue;

      if ((flags & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1)
      {
        representativeSubmatrixRow = submatrixRow;
        break;
      }
    }
    assert(representativeSubmatrixRow < SIZE_MAX);

    /* We then copy the nonzeros of that rank-1 subpart. */
    int8_t* representativeSubmatrixDense = NULL;
    size_t representativeSubmatrixColumn = SIZE_MAX;
    CMR_CALL( CMRallocStackArray(cmr, &representativeSubmatrixDense, sepa->numColumns) );
    for (size_t submatrixColumn = 0; submatrixColumn < sepa->numColumns; ++submatrixColumn)
      representativeSubmatrixDense[submatrixColumn] = 0;
    size_t first = matrix->rowSlice[submatrixRows[representativeSubmatrixRow]];
    size_t beyond = matrix->rowSlice[submatrixRows[representativeSubmatrixRow] + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t submatrixColumn = columnsToSubmatrixColumn[column];
      if (submatrixColumn == SIZE_MAX)
        continue;

      if ((sepa->columnsFlags[submatrixColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      {
        representativeSubmatrixDense[submatrixColumn] = matrix->entryValues[e];
        if (representativeSubmatrixColumn == SIZE_MAX)
          representativeSubmatrixColumn = submatrixColumn;
      }
    }

    /* We now go through all rows of that part and compare to the dense. */
    for (size_t submatrixRow = 0; submatrixRow < sepa->numRows; ++submatrixRow)
    {
      CMR_SEPA_FLAGS flags = sepa->rowsFlags[submatrixRow];

      /* Skip top-left rows, zero rows and the representative row. */
      if ((flags != (CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1)) || (submatrixRow == representativeSubmatrixRow))
        continue;

      size_t row = submatrixRows[submatrixRow];
      int8_t scaling = 0;
      first = matrix->rowSlice[row];
      beyond = matrix->rowSlice[row + 1];
      for (size_t e = first; e < beyond; ++e)
      {
        size_t column = matrix->entryColumns[e];
        size_t submatrixColumn = columnsToSubmatrixColumn[column];
        if (submatrixColumn == SIZE_MAX)
          continue;

        if ((sepa->columnsFlags[submatrixColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          continue;

        assert(representativeSubmatrixDense[submatrixColumn] != 0); /* Rank > 1 ??? */
        if (!scaling)
          scaling = representativeSubmatrixDense[submatrixColumn] * matrix->entryValues[e];

        if (matrix->entryValues[e] * scaling != representativeSubmatrixDense[submatrixColumn])
        {
          if (pviolator)
          {
            CMRdbgMsg(6, "-> not ternary!\n");
            CMR_CALL( CMRsubmatCreate(cmr, 2, 2, pviolator) );
            CMR_SUBMAT* violator = *pviolator;
            violator->rows[0] = representativeSubmatrixRow;
            violator->rows[1] = submatrixRow;
            violator->columns[0] = representativeSubmatrixColumn;
            violator->columns[1] = submatrixColumn;
          }
          *pisTernary = false;
          goto cleanup;
        }
      }
    }

cleanup:

    CMR_CALL( CMRfreeStackArray(cmr, &representativeSubmatrixDense) );
  }
  else
  {
    assert(false);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRsepaCheckTernary(CMR* cmr, CMR_SEPA* sepa, CMR_CHRMAT* matrix, bool* pisTernary,
  CMR_SUBMAT** pviolator)
{
  assert(cmr);
  assert(sepa);
  assert(matrix);
  assert(pisTernary);
  assert(!pviolator || !*pviolator);

  assert(matrix->numRows == sepa->numRows);
  assert(matrix->numColumns == sepa->numColumns);

  size_t* submatrixRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &submatrixRows, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    submatrixRows[row] = row;

  size_t* columnsSubmatrixColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnsSubmatrixColumn, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnsSubmatrixColumn[column] = column;

  CMR_CALL( checkTernary(cmr, sepa, matrix, submatrixRows, columnsSubmatrixColumn, pisTernary, pviolator) );

  CMR_CALL( CMRfreeStackArray(cmr, &columnsSubmatrixColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &submatrixRows) );

  return CMR_OKAY;
}

CMR_ERROR CMRsepaCheckTernarySubmatrix(CMR* cmr, CMR_SEPA* sepa, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix,
  bool* pisTernary, CMR_SUBMAT** pviolator)
{
  assert(cmr);
  assert(sepa);
  assert(matrix);
  assert(pisTernary);
  assert(!pviolator || !*pviolator);

  assert(sepa->numRows <= matrix->numRows);
  assert(sepa->numColumns <= matrix->numColumns);
  assert(sepa->numRows == submatrix->numRows);
  assert(sepa->numColumns == submatrix->numColumns);

  size_t* columnsSubmatrixColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnsSubmatrixColumn, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnsSubmatrixColumn[column] = SIZE_MAX;
  for (size_t submatrixColumn = 0; submatrixColumn < submatrix->numColumns; ++submatrixColumn)
    columnsSubmatrixColumn[submatrix->columns[submatrixColumn]] = submatrixColumn;

  CMR_CALL( checkTernary(cmr, sepa, matrix, submatrix->rows, columnsSubmatrixColumn, pisTernary, pviolator) );

  CMR_CALL( CMRfreeStackArray(cmr, &columnsSubmatrixColumn) );

  return CMR_OKAY;
}

CMR_ERROR CMRoneSumCompose(CMR* cmr, size_t numMatrices, CMR_CHRMAT** matrices, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(numMatrices > 0);
  assert(matrices);
  assert(presult);

  size_t numRows = 0;
  size_t numColumns = 0;
  size_t numNonzeros = 0;
  for (size_t i = 0; i < numMatrices; ++i)
  {
    numRows += matrices[i]->numRows;
    numColumns += matrices[i]->numColumns;
    numNonzeros += matrices[i]->numNonzeros;
  }

  CMR_CALL( CMRchrmatCreate(cmr, presult, numRows, numColumns, numNonzeros) );
  CMR_CHRMAT* result = *presult;

  size_t rows = 0;
  size_t columns = 0;
  size_t entry = 0;
  for (size_t i = 0; i < numMatrices; ++i)
  {
    CMR_CHRMAT* matrix = matrices[i];

    for (size_t r = 0; r < matrix->numRows; ++r)
    {
      result->rowSlice[rows + r] = entry;
      size_t beyond = matrix->rowSlice[r + 1];
      for (size_t e = matrix->rowSlice[r]; e < beyond; ++e)
      {
        result->entryColumns[entry] = columns + matrix->entryColumns[e];
        result->entryValues[entry++] = matrix->entryValues[e];
      }
    }

    rows += matrix->numRows;
    columns += matrix->numColumns;
  }

  result->rowSlice[rows] = entry;
  assert(rows == numRows);
  assert(columns == numColumns);
  assert(entry == numNonzeros);

  return CMR_OKAY;
}

CMR_ERROR CMRtwoSumCompose(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, size_t* firstSpecialRows,
  size_t* firstSpecialColumns, size_t* secondSpecialRows, size_t* secondSpecialColumns, int8_t characteristic,
  CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(first);
  assert(second);
  assert(presult);

  CMRdbgMsg(0, "CMRtwoSumCompose for a %zux%zu and a %zux%zu matrix.\n", first->numRows, first->numColumns,
    second->numRows, second->numColumns);

  size_t firstRowMarker = firstSpecialRows ? firstSpecialRows[0] : SIZE_MAX;
  size_t firstColumnMarker = firstSpecialColumns ? firstSpecialColumns[0] : SIZE_MAX;
  size_t secondRowMarker = secondSpecialRows ? secondSpecialRows[0] : SIZE_MAX;
  size_t secondColumnMarker = secondSpecialColumns ? secondSpecialColumns[0] : SIZE_MAX;
  bool bottomLeft;
  if ((firstRowMarker < SIZE_MAX) && (secondColumnMarker < SIZE_MAX)
    && (firstColumnMarker == SIZE_MAX) && (secondRowMarker == SIZE_MAX))
  {
    bottomLeft = true;
  }
  else if ((firstRowMarker == SIZE_MAX) && (secondColumnMarker == SIZE_MAX)
    && (firstColumnMarker < SIZE_MAX) && (secondRowMarker < SIZE_MAX))
  {
    bottomLeft = false;
  }
  else
    return CMR_ERROR_INPUT;

  char* markerColumn = NULL; /* Nonzero entries of the column vector among a,b. */
  size_t markerColumnNumNonzeros = 0; /* Number of nonzeros in markerColumn. */
  size_t markerRowNumNonzeros = 0; /* Number of nonzeros of a,b that is not markerColumn. */
  CMR_ERROR error = CMR_OKAY;

  /* Step 1: fill marker column. */
  if (bottomLeft)
  {
    /* rank 1 in bottom-right matrix. */
    markerRowNumNonzeros = first->rowSlice[firstRowMarker+1] - first->rowSlice[firstRowMarker];

    CMR_CALL( CMRallocStackArray(cmr, &markerColumn, second->numColumns) );
    for (size_t row = 0; row < second->numRows; ++row)
    {
      size_t entry;
      CMR_CALL( CMRchrmatFindEntry(second, row, secondColumnMarker, &entry) );
      if (entry == SIZE_MAX)
      {
        markerColumn[row] = 0;
        CMRdbgMsg(2, "markerColumn[%ld] is 0.\n", row);
      }
      else
      {
        assert(second->entryColumns[entry] == secondColumnMarker);
        markerColumn[row] = second->entryValues[entry];
        CMRdbgMsg(2, "markerColumn[%ld] = %d.\n", row, second->entryValues[entry]);
        ++markerColumnNumNonzeros;
      }
    }
  }
  else
  {
    /* rank 1 in top right */

    CMR_CALL( CMRallocStackArray(cmr, &markerColumn, first->numColumns) );
    for (size_t row = 0; row < first->numRows; ++row)
    {
      size_t entry;
      CMR_CALL( CMRchrmatFindEntry(first, row, firstColumnMarker, &entry) );
      if (entry == SIZE_MAX)
        markerColumn[row] = 0;
      else
      {
        assert(first->entryColumns[entry] == firstColumnMarker);
        markerColumn[row] = first->entryValues[entry];
        ++markerColumnNumNonzeros;
      }
    }

    markerRowNumNonzeros = second->rowSlice[secondRowMarker+1] - second->rowSlice[secondRowMarker];
  }

  /* Step 2: create the actual matrix. */
  CMRdbgMsg(0, "#nonzeros is %zu + %zu + %zu * %zu - 1\n", first->numNonzeros, second->numNonzeros,
    markerRowNumNonzeros - 1, markerColumnNumNonzeros - 1);
  CMR_CALL( CMRchrmatCreate(cmr, presult, first->numRows + second->numRows - 1,
    first->numColumns + second->numColumns - 1,
    first->numNonzeros + second->numNonzeros + (markerRowNumNonzeros - 1) * (markerColumnNumNonzeros - 1) - 1) );
  CMR_CHRMAT* result = *presult;
  size_t resultRow = 0;
  size_t resultNonzero = 0;

  /* Step 2a: top part. */
  for (size_t row = 0; row < first->numRows; ++row)
  {
    if (row != firstRowMarker)
    {
      result->rowSlice[resultRow] = resultNonzero;
      CMRdbgMsg(0, "Row %ld starts at nonzero #%ld.\n", resultRow, resultNonzero);

      /* First matrix is top-left. */
      for (size_t e = first->rowSlice[row]; e < first->rowSlice[row + 1]; ++e)
      {
        size_t column = first->entryColumns[e];
        if (column == firstColumnMarker)
          continue;
        result->entryValues[resultNonzero] = first->entryValues[e];
        result->entryColumns[resultNonzero] = column < firstColumnMarker ? column : (column - 1);
        CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
        ++resultNonzero;
      }

      /* Rank-1 matrix is top-right. */
      if (!bottomLeft && markerColumn[row])
      {
        for (size_t e = second->rowSlice[secondRowMarker]; e < second->rowSlice[secondRowMarker + 1]; ++e)
        {
          int resultValue = (int)(markerColumn[row]) * (int)second->entryValues[e];

          /* Take care of modulo operations. */
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }

          if (resultValue > INT8_MAX || resultValue < INT8_MIN)
          {
            error = CMR_ERROR_OVERFLOW;
            goto cleanup;
          }

          result->entryValues[resultNonzero] = resultValue;
          result->entryColumns[resultNonzero] = first->numColumns - 1 + second->entryColumns[e];
          CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
          ++resultNonzero;
        }
      }
      ++resultRow;
    }
  }

  /* Step 2b: bottom part. */
  for (size_t row = 0; row < second->numRows; ++row)
  {
    if (row != secondRowMarker)
    {
      result->rowSlice[resultRow] = resultNonzero;
      CMRdbgMsg(0, "Row %ld starts at nonzero #%ld.\n", resultRow, resultNonzero);

      /* Rank-1 matrix is bottom-left. */
      if (bottomLeft && markerColumn[row])
      {
        for (size_t e = first->rowSlice[firstRowMarker]; e < first->rowSlice[firstRowMarker + 1]; ++e)
        {
          size_t column = first->entryColumns[e];

          /* Take care of modulo operations. */
          int resultValue = (int)markerColumn[row] * (int)first->entryValues[e];
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }

          if (resultValue > INT8_MAX || resultValue < INT8_MIN)
          {
            error = CMR_ERROR_OVERFLOW;
            goto cleanup;
          }

          result->entryValues[resultNonzero] = resultValue;
          result->entryColumns[resultNonzero] = column;
          CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
          ++resultNonzero;
        }
      }
      
      /* Second matrix is bottom-right. */
      for (size_t e = second->rowSlice[row]; e < second->rowSlice[row + 1]; ++e)
      {
        size_t column = second->entryColumns[e];
        if (column == secondColumnMarker)
          continue;
        result->entryValues[resultNonzero] = second->entryValues[e];
        result->entryColumns[resultNonzero] = (bottomLeft ? first->numColumns : (first->numColumns - 1))
           + (column < secondColumnMarker ? column : (column - 1));
        CMRdbgMsg(2, "Added nonzero #%ld at %ld,%ld with value %d.\n", resultNonzero, resultRow,
          result->entryColumns[resultNonzero], result->entryValues[resultNonzero]);
        ++resultNonzero;
      }

      ++resultRow;
    }
  }
  result->rowSlice[result->numRows] = resultNonzero;
  assert(resultNonzero == result->numNonzeros);

  CMRdbgConsistencyAssert( CMRchrmatConsistency(result) );

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &markerColumn) );

  return error;
}

CMR_ERROR CMRtwoSumDecomposeFirst(CMR* cmr, CMR_CHRMAT* matrix, CMR_SEPA* sepa, CMR_CHRMAT** pfirst,
  size_t* firstRowsOrigin, size_t* firstColumnsOrigin, size_t* rowsToFirst, size_t* columnsToFirst,
  size_t* firstSpecialRows, size_t* firstSpecialColumns)
{
  assert(cmr);
  assert(matrix);
  assert(sepa);
  assert(sepa->type == CMR_SEPA_TYPE_TWO);
  assert(pfirst);

  CMRdbgMsg(0, "Computing first part of 2-sum decomposition of %zux%zu-matrix.\n", matrix->numRows, matrix->numColumns);

  /* Allocate missing arrays on stack. */
  bool hasFirstRowsOrigin = firstRowsOrigin;
  if (!hasFirstRowsOrigin)
    CMR_CALL( CMRallocStackArray(cmr, &firstRowsOrigin, matrix->numRows) );
  bool hasFirstColumnsOrigin = firstColumnsOrigin;
  if (!hasFirstColumnsOrigin)
    CMR_CALL( CMRallocStackArray(cmr, &firstColumnsOrigin, matrix->numColumns) );
  bool hasRowsToFirst = rowsToFirst;
  if (!hasRowsToFirst)
    CMR_CALL( CMRallocStackArray(cmr, &rowsToFirst, matrix->numRows) );
  bool hasColumnsToFirst = columnsToFirst;
  if (!hasColumnsToFirst)
    CMR_CALL( CMRallocStackArray(cmr, &columnsToFirst, matrix->numColumns) );

  char* denseColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseColumn, matrix->numColumns) );

  /* Number of rows of A. */
  size_t numRows = 0;
  size_t extraRow = SIZE_MAX;
  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    if ((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
    {
      rowsToFirst[row] = numRows;
      firstRowsOrigin[numRows++] = row;
    }
    else
    {
      rowsToFirst[row] = SIZE_MAX;
      if (extraRow == SIZE_MAX && sepa->rowsFlags[row] & CMR_SEPA_MASK_EXTRA)
        extraRow = row;
    }
  }

  /* Number of columns of A. */
  size_t numColumns = 0;
  size_t extraColumn = SIZE_MAX;
  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    if ((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
    {
      columnsToFirst[column] = numColumns;
      firstColumnsOrigin[numColumns++] = column;
    }
    else
    {
      columnsToFirst[column] = SIZE_MAX;
      if (extraColumn == SIZE_MAX && sepa->columnsFlags[column] & CMR_SEPA_MASK_EXTRA)
        extraColumn = column;
    }
  }

  if (extraRow < SIZE_MAX)
    firstRowsOrigin[numRows++] = extraRow;
  else
    firstColumnsOrigin[numColumns++] = extraColumn;

  /* Count number of nonzeros and copy column vector. */
  size_t numNonzeros = 0;
  for (size_t row1 = 0; row1 < numRows; ++row1)
  {
    size_t row = firstRowsOrigin[row1];
    if (row == SIZE_MAX)
      row = extraRow;

    if (extraColumn < SIZE_MAX)
      denseColumn[row1] = 0;

    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t column1 = columnsToFirst[column];
      if (column1 < SIZE_MAX)
        numNonzeros++;
      else if (column == extraColumn)
      {
        numNonzeros++;
        denseColumn[row1] = matrix->entryValues[e];
      }
    }
  }

  /* Copy the matrix entries. */
  CMR_CALL( CMRchrmatCreate(cmr, pfirst, numRows, numColumns, numNonzeros) );
  CMR_CHRMAT* first = *pfirst;
  first->numNonzeros = 0;
  for (size_t row1 = 0; row1 < numRows; ++row1)
  {
    size_t row = firstRowsOrigin[row1];
    if (row == SIZE_MAX)
      row = extraRow;

    /* Row from A and potentially c^T. */
    first->rowSlice[row1] = first->numNonzeros;
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t column1 = columnsToFirst[column];
      if (column1 < SIZE_MAX)
      {
        first->entryColumns[first->numNonzeros] = column1;
        first->entryValues[first->numNonzeros++] = matrix->entryValues[e];
      }
    }

    /* Nonzero of a. */
    if (extraColumn < SIZE_MAX && denseColumn[row1])
    {
      first->entryColumns[first->numNonzeros] = numColumns - 1;
      first->entryValues[first->numNonzeros++] = denseColumn[row1];
    }
  }
  first->rowSlice[first->numRows] = first->numNonzeros;

  /* Set marker to last row or column. */
  if (firstSpecialRows)
    firstSpecialRows[0] = (extraRow < SIZE_MAX) ? (numRows - 1) : SIZE_MAX;
  if (firstSpecialColumns)
    firstSpecialColumns[0] = (extraRow < SIZE_MAX) ? SIZE_MAX : (numColumns - 1);

  /* Free local arrays. */
  CMR_CALL( CMRfreeStackArray(cmr, &denseColumn) );
  if (hasColumnsToFirst)
    CMR_CALL( CMRfreeStackArray(cmr, &columnsToFirst) );
  if (hasRowsToFirst)
    CMR_CALL( CMRfreeStackArray(cmr, &rowsToFirst) );
  if (hasFirstColumnsOrigin)
    CMR_CALL( CMRfreeStackArray(cmr, &firstColumnsOrigin) );
  if (hasFirstRowsOrigin)
    CMR_CALL( CMRfreeStackArray(cmr, &firstRowsOrigin) );

  return CMR_OKAY;
}


CMR_ERROR CMRtwoSumDecomposeSecond(CMR* cmr, CMR_CHRMAT* matrix, CMR_SEPA* sepa, CMR_CHRMAT** psecond,
  size_t* secondRowsOrigin, size_t* secondColumnsOrigin, size_t* rowsToSecond, size_t* columnsToSecond,
  size_t* secondSpecialRows, size_t* secondSpecialColumns)
{
  assert(cmr);
  assert(matrix);
  assert(sepa);
  assert(sepa->type == CMR_SEPA_TYPE_TWO);
  assert(psecond);

  CMRdbgMsg(0, "Computing second part of 2-sum decomposition of %zux%zu-matrix.\n", matrix->numRows, matrix->numColumns);

  /* Allocate missing arrays on stack. */
  bool hasSecondRowsOrigin = secondRowsOrigin;
  if (!hasSecondRowsOrigin)
    CMR_CALL( CMRallocStackArray(cmr, &secondRowsOrigin, matrix->numRows) );
  bool hasSecondColumnsOrigin = secondColumnsOrigin;
  if (!hasSecondColumnsOrigin)
    CMR_CALL( CMRallocStackArray(cmr, &secondColumnsOrigin, matrix->numColumns) );
  bool hasRowsToSecond = rowsToSecond;
  if (!hasRowsToSecond)
    CMR_CALL( CMRallocStackArray(cmr, &rowsToSecond, matrix->numRows) );
  bool hasColumnsToSecond = columnsToSecond;
  if (!hasColumnsToSecond)
    CMR_CALL( CMRallocStackArray(cmr, &columnsToSecond, matrix->numColumns) );

  char* denseColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseColumn, matrix->numColumns) );

  /* Find extra row. */
  size_t extraRow = SIZE_MAX;
  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    if ((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST && (sepa->rowsFlags[row] & CMR_SEPA_MASK_EXTRA))
    {
      extraRow = row;
      break;
    }
  }
  size_t extraColumn = SIZE_MAX;
  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    if ((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST
      && (sepa->columnsFlags[column] & CMR_SEPA_MASK_EXTRA))
    {
      extraColumn = column;
      break;
    }
  }

  CMRdbgMsg(2, "Extra row = %zu, extra column = %zu\n", extraRow, extraColumn);

  /* Number of rows of D. */
  size_t numRows = extraRow < SIZE_MAX ? 1 : 0;
  secondRowsOrigin[0] = extraRow;
  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    if ((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
    {
      rowsToSecond[row] = numRows;
      secondRowsOrigin[numRows++] = row;
    }
    else
      rowsToSecond[row] = SIZE_MAX;
  }

  /* Number of columns of D. */
  size_t numColumns = extraColumn < SIZE_MAX ? 1 : 0;
  secondColumnsOrigin[0] = extraColumn;
  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    if ((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
    {
      columnsToSecond[column] = numColumns;
      secondColumnsOrigin[numColumns++] = column;
    }
    else
      columnsToSecond[column] = SIZE_MAX;
  }

  /* Count number of nonzeros, copy column vector, and find out if we need to negate it. */
  size_t numNonzeros = 0;
  char scale = 0;

  for (size_t row2 = 0; row2 < numRows; ++row2)
  {
    size_t row = (extraRow < SIZE_MAX && row2 == 0) ? extraRow : secondRowsOrigin[row2];

    if (extraColumn < SIZE_MAX)
      denseColumn[row2] = 0;

    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t column2 = columnsToSecond[column];
      if (column2 < SIZE_MAX)
      {
        numNonzeros++;
        if (row == extraRow && scale == 0)
        {
          CMRdbgMsg(2, "Scaling entry is %d, taken from row r%zu, c%zu.\n", matrix->entryValues[e], row+1, column+1);
          scale = matrix->entryValues[e];
        }
      }
      else if (column == extraColumn)
      {
        numNonzeros++;
        denseColumn[row2] = matrix->entryValues[e];
        if (scale == 0)
        {
          scale = matrix->entryValues[e];
          CMRdbgMsg(2, "Scaling entry is %d, taken from row r%zu, c%zu.\n", matrix->entryValues[e], row+1, column+1);
        }
      }
    }
  }
  assert(scale != 0);

  CMRdbgMsg(2, "Scaling rank-1 vector by %d.\n", scale);

  /* Copy the matrix entries. */
  CMR_CALL( CMRchrmatCreate(cmr, psecond, numRows, numColumns, numNonzeros) );
  CMR_CHRMAT* second = *psecond;
  second->numNonzeros = 0;
  for (size_t row2 = 0; row2 < numRows; ++row2)
  {
    size_t row = (extraRow < SIZE_MAX && row2 == 0) ? extraRow : secondRowsOrigin[row2];
    second->rowSlice[row2] = second->numNonzeros;

    /* Nonzero from d. */
    if (extraColumn < SIZE_MAX && denseColumn[row2])
    {
      second->entryColumns[second->numNonzeros] = 0;
      second->entryValues[second->numNonzeros++] = scale * denseColumn[row2];
    }

    /* Row from [ b^T \\ D ]. */
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t column2 = columnsToSecond[column];
      if (column2 < SIZE_MAX)
      {
        char value = matrix->entryValues[e];
        if (row == extraRow)
          value *= scale;
        second->entryColumns[second->numNonzeros] = column2;
        second->entryValues[second->numNonzeros++] = value;
      }
    }
  }
  second->rowSlice[second->numRows] = second->numNonzeros;

  /* Set marker to last row or column. */
  if (secondSpecialRows)
    secondSpecialRows[0] = (extraRow < SIZE_MAX) ? 0 : SIZE_MAX;
  if (secondSpecialColumns)
    secondSpecialColumns[0] = (extraRow < SIZE_MAX) ? SIZE_MAX : 0;

  /* Free local arrays. */
  CMR_CALL( CMRfreeStackArray(cmr, &denseColumn) );
  if (hasColumnsToSecond)
    CMR_CALL( CMRfreeStackArray(cmr, &columnsToSecond) );
  if (hasRowsToSecond)
    CMR_CALL( CMRfreeStackArray(cmr, &rowsToSecond) );
  if (hasSecondColumnsOrigin)
    CMR_CALL( CMRfreeStackArray(cmr, &secondColumnsOrigin) );
  if (hasSecondRowsOrigin)
    CMR_CALL( CMRfreeStackArray(cmr, &secondRowsOrigin) );

  return CMR_OKAY;
}


CMR_ERROR CMRthreeSumSeymourCompose(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, size_t* firstSpecialRows,
  size_t* firstSpecialColumns, size_t* secondSpecialRows, size_t* secondSpecialColumns, int8_t characteristic,
  CMR_CHRMAT ** presult)
{
  assert(cmr);
  assert(first);
  assert(second);
  assert(presult);

  CMRdbgMsg(0, "CMRthreeSumSeymourCompose for a %zux%zu and a %zux%zu matrix.\n", first->numRows, first->numColumns,
    second->numRows, second->numColumns);

  if (!firstSpecialRows || (firstSpecialRows[0] >= first->numRows))
    return CMR_ERROR_INPUT;
  if (!firstSpecialColumns || (firstSpecialColumns[0] >= first->numColumns)
    || (firstSpecialColumns[1] >= first->numColumns))
  {
    return CMR_ERROR_INPUT;
  }
  if (!secondSpecialRows || (secondSpecialRows[0] >= second->numRows))
    return CMR_ERROR_INPUT;
  if (!secondSpecialColumns || (secondSpecialColumns[0] >= second->numColumns)
    || (secondSpecialColumns[1] >= second->numColumns))
  {
    return CMR_ERROR_INPUT;
  }

  size_t firstSpecialColumnLeft;
  size_t firstSpecialColumnRight;
  if (firstSpecialColumns[0] < firstSpecialColumns[1])
  {
    firstSpecialColumnLeft = firstSpecialColumns[0];
    firstSpecialColumnRight = firstSpecialColumns[1];
  }
  else if (firstSpecialColumns[0] > firstSpecialColumns[1])
  {
    firstSpecialColumnLeft = firstSpecialColumns[1];
    firstSpecialColumnRight = firstSpecialColumns[0];
  }
  else
    return CMR_ERROR_INPUT;

  size_t secondSpecialColumnLeft;
  size_t secondSpecialColumnRight;
  if (secondSpecialColumns[0] < secondSpecialColumns[1])
  {
    secondSpecialColumnLeft = secondSpecialColumns[0];
    secondSpecialColumnRight = secondSpecialColumns[1];
  }
  else if (secondSpecialColumns[0] > secondSpecialColumns[1])
  {
    secondSpecialColumnLeft = secondSpecialColumns[1];
    secondSpecialColumnRight = secondSpecialColumns[0];
  }
  else
    return CMR_ERROR_INPUT;


  /* Number of nonzeros. */
  size_t firstMainNumNonzeros = 0;
  size_t secondMainNumNonzeros = 0;
  size_t firstSpecialRowNumNonzeros = first->rowSlice[firstSpecialRows[0] + 1] - first->rowSlice[firstSpecialColumns[0]] - 1;
  size_t secondSpecialRowNumNonzeros = second->rowSlice[secondSpecialRows[0] + 1] - second->rowSlice[secondSpecialRows[0]] - 1;
  char* firstSpecialColumnDense = NULL; /* Entries of a. */
  size_t firstSpecialColumnNumNonzeros = 0;
  char* secondSpecialColumnDense = NULL; /* Entries of d. */
  size_t secondSpecialColumnNumNonzeros = 0;
  size_t specialColumn1CopyNumNonzeros; /* Number of nonzeros in copy of a or copy of d. */
  CMR_ERROR error = CMR_OKAY;
  CMR_CALL( CMRallocStackArray(cmr, &firstSpecialColumnDense, first->numRows - 1) );
  CMR_CALL( CMRallocStackArray(cmr, &secondSpecialColumnDense, second->numRows - 1) );
  size_t firstMainRow = 0;
  specialColumn1CopyNumNonzeros = 0;
  char firstEpsilon = 0;
  for (size_t firstRow = 0; firstRow < first->numRows; ++firstRow)
  {
    if (firstRow != firstSpecialRows[0])
      firstSpecialColumnDense[firstMainRow] = 0;
    size_t beyond = first->rowSlice[firstRow + 1];
    for (size_t e = first->rowSlice[firstRow]; e < beyond; ++e)
    {
      size_t firstColumn = first->entryColumns[e];
      if (firstRow != firstSpecialRows[0])
      {
        if (firstColumn == firstSpecialColumnLeft)
        {
          firstSpecialColumnDense[firstMainRow] = first->entryValues[e];
          ++firstSpecialColumnNumNonzeros;
        }
        else if (firstColumn == firstSpecialColumnRight)
        {
          if (firstSpecialColumnDense[firstMainRow] != first->entryValues[e])
          {
            CMRdbgMsg(2, "Special columns of 1st matrix differ: r%zu,c%zu -> %d, r%zu,c%zu -> %d.\n",
              firstRow+1, firstSpecialColumnLeft+1, firstSpecialColumnDense[firstMainRow], firstRow+1,
              firstSpecialColumnRight+1, first->entryValues[e]);
            error = CMR_ERROR_STRUCTURE;
            goto cleanup;
          }
          ++specialColumn1CopyNumNonzeros;
        }
        else
        {
          ++firstMainNumNonzeros;
        }
      }
      else
      {
        if (firstColumn == firstSpecialColumns[0])
        {
          CMRdbgMsg(2, "Expected a 0-entry in 1st matrix: r%zu,c%zu, but %d found.\n",
              firstRow+1, firstSpecialColumns[0]+1, first->entryValues[e]);
            error = CMR_ERROR_STRUCTURE;
            goto cleanup;
        }
        else if (firstColumn == firstSpecialColumns[1])
          firstEpsilon = first->entryValues[e];
      }
    }

    if (firstRow != firstSpecialRows[0])
      ++firstMainRow;
  }

  /* Check special columns of 1st matrix. */
  if (specialColumn1CopyNumNonzeros != firstSpecialColumnNumNonzeros)
  {
    CMRdbgMsg(2, "Number of nonzeros in special columns of 1st matrix is %zu and %zu.\n", firstSpecialColumnNumNonzeros,
      specialColumn1CopyNumNonzeros);
    error = CMR_ERROR_STRUCTURE;
    goto cleanup;
  }
  else if (firstEpsilon == 0)
  {
    CMRdbgMsg(2, "Epsilon-entry of 1st matrix is %d.\n", firstEpsilon);
    error = CMR_ERROR_STRUCTURE;
    goto cleanup;
  }

  specialColumn1CopyNumNonzeros = 0;
  size_t secondMainRow = 0;
  char secondEpsilon = 0;
  for (size_t secondRow = 0; secondRow < second->numRows; ++secondRow)
  {
    if (secondRow != secondSpecialRows[0])
      secondSpecialColumnDense[secondMainRow] = 0;
    size_t beyond = second->rowSlice[secondRow + 1];
    for (size_t e = second->rowSlice[secondRow]; e < beyond; ++e)
    {
      size_t secondColumn = second->entryColumns[e];
      if (secondRow != secondSpecialRows[0])
      {
        if (secondColumn == secondSpecialColumnLeft)
        {
          secondSpecialColumnDense[secondMainRow] = second->entryValues[e];
          ++secondSpecialColumnNumNonzeros;
        }
        else if (secondColumn == secondSpecialColumnRight)
        {
          if (secondSpecialColumnDense[secondMainRow] != second->entryValues[e])
          {
            CMRdbgMsg(2, "Special columns of 2nd matrix differ: r%zu,c%zu -> %d, r%zu,c%zu -> %d.\n",
              secondRow+1, secondSpecialColumnLeft+1, secondSpecialColumnDense[secondMainRow], secondRow+1,
              secondSpecialColumnRight+1, second->entryValues[e]);
            error = CMR_ERROR_STRUCTURE;
            goto cleanup;
          }
          ++specialColumn1CopyNumNonzeros;
        }
        else
        {
          ++secondMainNumNonzeros;
        }
      }
      else
      {
        if (secondColumn == secondSpecialColumns[1])
        {
          CMRdbgMsg(2, "Expected a 0-entry in 2nd matrix: r%zu,c%zu, but %d found.\n",
              secondRow+1, secondSpecialColumns[1]+1, second->entryValues[e]);
            error = CMR_ERROR_STRUCTURE;
            goto cleanup;
        }
        else if (secondColumn == secondSpecialColumns[0])
          secondEpsilon = second->entryValues[e];
      }
    }

    if (secondRow != secondSpecialRows[0])
      ++secondMainRow;
  }

  /* Check special columns. */
  if (specialColumn1CopyNumNonzeros != secondSpecialColumnNumNonzeros)
  {
    CMRdbgMsg(2, "Number of nonzeros in special columns of 2nd matrix is %zu and %zu.\n",
      secondSpecialColumnNumNonzeros, specialColumn1CopyNumNonzeros);
    error = CMR_ERROR_STRUCTURE;
    goto cleanup;
  }
  else if (secondEpsilon == 0)
  {
    CMRdbgMsg(2, "Epsilon-entry of 2nd matrix is 0.\n");
    error = CMR_ERROR_STRUCTURE;
    goto cleanup;
  }
  else if (secondEpsilon != firstEpsilon)
  {
    CMRdbgMsg(2, "Epsilon-entries of matrices are %d and %d.\n", firstEpsilon, secondEpsilon);
    error = CMR_ERROR_STRUCTURE;
    goto cleanup;
  }

  CMRdbgMsg(2, "First: main has %zu nonzeros, column has %zu nonzeros and row has %zu nonzeros.\n",
    firstMainNumNonzeros, firstSpecialColumnNumNonzeros, firstSpecialRowNumNonzeros);
  CMRdbgMsg(2, "Second: main has %zu nonzeros, column has %zu nonzeros and row has %zu nonzeros.\n",
    secondMainNumNonzeros, secondSpecialColumnNumNonzeros, secondSpecialRowNumNonzeros);

  /* Create resulting matrix. */
  CMR_CALL( CMRchrmatCreate(cmr, presult, first->numRows + second->numRows - 2,
    first->numColumns + second->numColumns - 4, firstMainNumNonzeros + secondMainNumNonzeros
    + firstSpecialRowNumNonzeros * secondSpecialColumnNumNonzeros
    + secondSpecialRowNumNonzeros * firstSpecialColumnNumNonzeros) );
  CMR_CHRMAT* result = *presult;

  CMRdbgMsg(2, "Resulting matrix is %zux%zu with %zu nonzeros.\n", result->numRows, result->numColumns,
    result->numNonzeros);

  size_t row = 0;
  size_t entry = 0;
  for (size_t firstRow = 0; firstRow < first->numRows; ++firstRow)
  {
    result->rowSlice[row] = entry;
    if (firstRow != firstSpecialRows[0])
    {
      /* Nonzeros in top-left. */
      size_t beyond = first->rowSlice[firstRow+1];
      for (size_t e = first->rowSlice[firstRow]; e < beyond; ++e)
      {
        size_t column = first->entryColumns[e];
        if (column == firstSpecialColumnLeft || column == firstSpecialColumnRight)
          continue;
        else if (column > firstSpecialColumnRight)
          column -= 2;
        else if (column > firstSpecialColumnLeft)
          column -= 1;

        CMRdbgMsg(4, "First main nonzero r%zu,c%zu -> %d is mapped to r%zu,c%zu.\n", firstRow+1,
          first->entryColumns[e]+1, first->entryValues[e], row+1, column+1);
        result->entryColumns[entry] = column;
        result->entryValues[entry] = first->entryValues[e];
        ++entry;
      }

      /* Nonzeros in top-right. */
      if (firstSpecialColumnDense[row])
      {
        char factor = firstSpecialColumnDense[row];
        beyond = second->rowSlice[secondSpecialRows[0] + 1];
        for (size_t e = second->rowSlice[secondSpecialRows[0]]; e < beyond; ++e)
        {
          size_t column = second->entryColumns[e];
          if (column == secondSpecialColumnLeft || column == secondSpecialColumnRight)
            continue;
          else if (column > secondSpecialColumnRight)
            column -= 2;
          else if (column > secondSpecialColumnLeft)
            column -= 1;

          CMRdbgMsg(4, "Top-right nonzero r%zu (first) and c%zu (second) -> %d*%d is mapped to r%zu,c%zu.\n", firstRow+1,
            second->entryColumns[e]+1, factor, second->entryValues[e], row+1, first->numColumns - 2 + column+1);
          result->entryColumns[entry] = first->numColumns - 2 + column;
          int resultValue = (int) factor * (int) second->entryValues[e];
          assert(resultValue != 0);
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }
          assert(resultValue != 0);
          result->entryValues[entry] = resultValue;
          ++entry;
        }
      }
      ++row;
    }
  }
  secondMainRow = 0;
  for (size_t secondRow = 0; secondRow < second->numRows; ++secondRow)
  {
    result->rowSlice[row] = entry;
    if (secondRow != secondSpecialRows[0])
    {
      /* Nonzeros in bottom-left. */
      if (secondSpecialColumnDense[secondMainRow])
      {
        char factor = secondSpecialColumnDense[secondMainRow];
        size_t beyond = first->rowSlice[firstSpecialRows[0] + 1];
        for (size_t e = first->rowSlice[firstSpecialRows[0]]; e < beyond; ++e)
        {
          size_t column = first->entryColumns[e];
          if (column == firstSpecialColumnLeft || column == firstSpecialColumnRight)
            continue;
          else if (column > firstSpecialColumnRight)
            column -= 2;
          else if (column > firstSpecialColumnLeft)
            column -= 1;

          CMRdbgMsg(4, "Bottom-left nonzero r%zu (second) and c%zu (first) -> %d*%d is mapped to r%zu,c%zu.\n",
            secondRow+1, first->entryColumns[e]+1, factor, first->entryValues[e], row+1, column+1);
          result->entryColumns[entry] = column;
          int resultValue = (int) factor * (int) first->entryValues[e];
          assert(resultValue != 0);
          if (characteristic != 0)
          {
            resultValue = resultValue % characteristic;
            if (resultValue < 0)
              resultValue += characteristic;
            if (characteristic == 3 && resultValue == 2)
              resultValue -= 3;
          }
          assert(resultValue != 0);
          result->entryValues[entry] = resultValue;
          ++entry;
        }
      }

      /* Nonzeros in bottom-right. */
      size_t beyond = second->rowSlice[secondRow+1];
      for (size_t e = second->rowSlice[secondRow]; e < beyond; ++e)
      {
        size_t column = second->entryColumns[e];
        if (column == secondSpecialColumnLeft || column == secondSpecialColumnRight)
          continue;
        else if (column > secondSpecialColumnRight)
          column -= 2;
        else if (column > secondSpecialColumnLeft)
          column -= 1;

        CMRdbgMsg(4, "Second main nonzero r%zu,c%zu -> %d is mapped to r%zu,c%zu.\n", secondRow+1,
          second->entryColumns[e]+1, second->entryValues[e], row+1, first->numColumns - 2 + column+1);
        result->entryColumns[entry] = first->numColumns - 2 + column;
        result->entryValues[entry] = second->entryValues[e];
        ++entry;
      }

      ++row;
      ++secondMainRow;
    }
  }
  result->rowSlice[row] = entry;
  assert(entry == result->numNonzeros);

cleanup:
  CMR_CALL( CMRfreeStackArray(cmr, &secondSpecialColumnDense) );
  CMR_CALL( CMRfreeStackArray(cmr, &firstSpecialColumnDense) );

  return error;
}

CMR_ERROR CMRthreeSumSeymourDecomposeEpsilon(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose, CMR_SEPA* sepa,
  char* pepsilon)
{
  assert(cmr);
  assert(matrix);
  assert(transpose);
  assert(sepa);
  assert(pepsilon);

  size_t rightRow = SIZE_MAX;
  size_t bottomColumn = SIZE_MAX;

  /* Find first nonzero row in top-right part and first nonzero column in bottom-left part. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      if (rightRow == SIZE_MAX && (sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST
        && (sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
      {
        rightRow = row;
      }
      if (bottomColumn == SIZE_MAX && (sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND
        && (sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      {
        bottomColumn = column;
      }
    }
  }

  /* Graph search in top-left matrix: rows/columns of second are disabled (group -1), rank-1 rows are sources (group 1),
   * rank-1 columns are targets (group 2), all other rows/columns of first are neutral (group 0). */

  bool connected;
  CMR_ELEMENT source;
  CMR_ELEMENT target;
  int sum;
  int* rowsGroup = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowsGroup, matrix->numRows) );
  int* columnsGroup = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnsGroup, matrix->numColumns) );
  size_t rowBottomLeft = SIZE_MAX;
  size_t columnTopRight = SIZE_MAX;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    rowsGroup[row] = ((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND) ? -1 :
      ( (sepa->rowsFlags[row] & CMR_SEPA_MASK_EXTRA) ? 1 : 0 );
    if (rowBottomLeft == SIZE_MAX && ((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
      && ((sepa->rowsFlags[row] & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1))
      rowBottomLeft = row;
  }
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    columnsGroup[column] = ((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND) ? -1 :
      ( (sepa->columnsFlags[column] & CMR_SEPA_MASK_EXTRA) ? 2 : 0 );
    if (columnTopRight == SIZE_MAX && ((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
      && ((sepa->columnsFlags[column] & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1))
      columnTopRight = column;
  }

  if (rowBottomLeft == SIZE_MAX || columnTopRight == SIZE_MAX)
  {
    CMR_CALL( CMRfreeStackArray(cmr, &columnsGroup) );
    CMR_CALL( CMRfreeStackArray(cmr, &rowsGroup) );
    return CMR_ERROR_INPUT;
  }

  CMR_CALL( CMRchrmatSubmatrixBipartitePath(cmr, matrix, transpose, rowsGroup, columnsGroup, &connected, &source,
    &target, 0, 0, &sum) );

  CMR_CALL( CMRfreeStackArray(cmr, &columnsGroup) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowsGroup) );

  if (!connected)
    return CMR_ERROR_INPUT;

  /* Find entry of source row in first top-left rank-1 column. */
  assert(CMRelementIsRow(source));
  size_t sourceRow = CMRelementToRowIndex(source);
  size_t beyond = matrix->rowSlice[sourceRow + 1];
  char entryTopRight = 0;
  for (size_t e = matrix->rowSlice[sourceRow]; e < beyond; ++e)
  {
    if (matrix->entryColumns[e] == columnTopRight)
    {
      entryTopRight = matrix->entryValues[e];
      break;
    }
  }

  /* Find entry of target column in first bottom-right rank-1 row. */
  assert(CMRelementIsColumn(target));
  beyond = matrix->rowSlice[rowBottomLeft + 1];
  char entryBottomLeft = 0;
  for (size_t e = matrix->rowSlice[rowBottomLeft]; e < beyond; ++e)
  {
    if (matrix->entryColumns[e] == CMRelementToColumnIndex(target))
    {
      entryBottomLeft = matrix->entryValues[e];
      break;
    }
  }
  assert(entryBottomLeft);

  /* sum + entries in a and c^T. */
  *pepsilon = ((sum + entryTopRight + entryBottomLeft) % 4 == 1) ? -1 : +1;

  assert((sum + entryTopRight + entryBottomLeft + *pepsilon) % 4 == 0);

  return CMR_OKAY;
}


CMR_ERROR CMRthreeSumSeymourDecomposeFirst(CMR* cmr, CMR_CHRMAT* matrix, CMR_SEPA* sepa, char epsilon,
  CMR_CHRMAT** pfirst, size_t* firstRowsOrigin, size_t* firstColumnsOrigin, size_t* rowsToFirst, size_t* columnsToFirst,
  size_t* firstSpecialRows, size_t* firstSpecialColumns)
{
  assert(cmr);
  assert(matrix);
  assert(sepa);
  assert(epsilon == -1 || epsilon == +1);
  assert(pfirst);

  CMRdbgMsg(0, "Computing first part of Seymour 3-sum decomposition of %zux%zu-matrix.\n", matrix->numRows,
    matrix->numColumns);

  /* Allocate missing arrays on stack. */
  bool hasFirstRowsOrigin = firstRowsOrigin;
  if (!hasFirstRowsOrigin)
    CMR_CALL( CMRallocStackArray(cmr, &firstRowsOrigin, matrix->numRows) );
  bool hasFirstColumnsOrigin = firstColumnsOrigin;
  if (!hasFirstColumnsOrigin)
    CMR_CALL( CMRallocStackArray(cmr, &firstColumnsOrigin, matrix->numColumns) );
  bool hasRowsToFirst = rowsToFirst;
  if (!hasRowsToFirst)
    CMR_CALL( CMRallocStackArray(cmr, &rowsToFirst, matrix->numRows) );
  bool hasColumnsToFirst = columnsToFirst;
  if (!hasColumnsToFirst)
    CMR_CALL( CMRallocStackArray(cmr, &columnsToFirst, matrix->numColumns) );

  char* denseColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseColumn, matrix->numColumns) );

  /* Number of rows of A. */
  size_t numRows = 0;
  size_t extraRow = SIZE_MAX;
  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    if ((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
    {
      rowsToFirst[row] = numRows;
      firstRowsOrigin[numRows++] = row;
    }
    else
    {
      rowsToFirst[row] = SIZE_MAX;
      if (extraRow == SIZE_MAX && sepa->rowsFlags[row] & CMR_SEPA_MASK_EXTRA)
        extraRow = row;
    }
  }

  /* Number of columns of A. */
  size_t numColumns = 0;
  size_t extraColumn = SIZE_MAX;
  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    if ((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
    {
      columnsToFirst[column] = numColumns;
      firstColumnsOrigin[numColumns++] = column;
    }
    else
    {
      columnsToFirst[column] = SIZE_MAX;
      if (extraColumn == SIZE_MAX && sepa->columnsFlags[column] & CMR_SEPA_MASK_EXTRA)
        extraColumn = column;
    }
  }

  assert(extraRow < SIZE_MAX);
  assert(extraColumn < SIZE_MAX);

  firstRowsOrigin[numRows++] = extraRow;
  firstColumnsOrigin[numColumns++] = extraColumn;
  firstColumnsOrigin[numColumns++] = extraColumn;

  /* Count number of nonzeros and copy column vector. */
  size_t numNonzeros = 1;
  for (size_t row1 = 0; row1 < numRows; ++row1)
  {
    size_t row = firstRowsOrigin[row1];
    denseColumn[row1] = 0;
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t column1 = columnsToFirst[column];
      if (column1 < SIZE_MAX)
        numNonzeros++;
      else if (column == extraColumn && row != extraRow)
      {
        numNonzeros += 2;
        denseColumn[row1] = matrix->entryValues[e];
      }
    }
  }

  /* Copy the matrix entries. */
  CMR_CALL( CMRchrmatCreate(cmr, pfirst, numRows, numColumns, numNonzeros) );
  CMR_CHRMAT* first = *pfirst;
  first->numNonzeros = 0;
  for (size_t row1 = 0; row1 < numRows; ++row1)
  {
    size_t row = firstRowsOrigin[row1];

    /* Row from A and potentially c^T. */
    first->rowSlice[row1] = first->numNonzeros;
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t column1 = columnsToFirst[column];
      if (column1 < SIZE_MAX)
      {
        first->entryColumns[first->numNonzeros] = column1;
        first->entryValues[first->numNonzeros++] = matrix->entryValues[e];
      }
    }

    if (row != extraRow)
    {
      /* Nonzero of a. */
      if (denseColumn[row1])
      {
        first->entryColumns[first->numNonzeros] = numColumns - 2;
        first->entryValues[first->numNonzeros++] = denseColumn[row1];
        first->entryColumns[first->numNonzeros] = numColumns - 1;
        first->entryValues[first->numNonzeros++] = denseColumn[row1];
      }
    }
    else
    {
      /* Nonzero below a. */
      first->entryColumns[first->numNonzeros] = numColumns - 1;
      first->entryValues[first->numNonzeros++] = epsilon;
    }
  }
  first->rowSlice[first->numRows] = first->numNonzeros;

  /* Set markers to last row and last two columns. */
  if (firstSpecialRows)
    firstSpecialRows[0] = numRows - 1;
  if (firstSpecialColumns)
  {
    firstSpecialColumns[0] = numColumns - 2;
    firstSpecialColumns[1] = numColumns - 1;
  }

  /* Free local arrays. */
  CMR_CALL( CMRfreeStackArray(cmr, &denseColumn) );
  if (hasColumnsToFirst)
    CMR_CALL( CMRfreeStackArray(cmr, &columnsToFirst) );
  if (hasRowsToFirst)
    CMR_CALL( CMRfreeStackArray(cmr, &rowsToFirst) );
  if (hasFirstColumnsOrigin)
    CMR_CALL( CMRfreeStackArray(cmr, &firstColumnsOrigin) );
  if (hasFirstRowsOrigin)
    CMR_CALL( CMRfreeStackArray(cmr, &firstRowsOrigin) );

  return CMR_OKAY;
}

CMR_ERROR CMRthreeSumSeymourDecomposeSecond(CMR* cmr, CMR_CHRMAT* matrix, CMR_SEPA* sepa, char epsilon,
  CMR_CHRMAT** psecond, size_t* secondRowsOrigin, size_t* secondColumnsOrigin, size_t* rowsToSecond,
  size_t* columnsToSecond, size_t* secondSpecialRows, size_t* secondSpecialColumns)
{
  assert(cmr);
  assert(matrix);
  assert(sepa);
  assert(epsilon == -1 || epsilon == +1);
  assert(psecond);

  CMRdbgMsg(0, "Computing second part of Seymour 3-sum decomposition of %zux%zu-matrix.\n", matrix->numRows,
    matrix->numColumns);

  /* Allocate missing arrays on stack. */
  bool hasSecondRowsOrigin = secondRowsOrigin;
  if (!hasSecondRowsOrigin)
    CMR_CALL( CMRallocStackArray(cmr, &secondRowsOrigin, matrix->numRows) );
  bool hasSecondColumnsOrigin = secondColumnsOrigin;
  if (!hasSecondColumnsOrigin)
    CMR_CALL( CMRallocStackArray(cmr, &secondColumnsOrigin, matrix->numColumns) );
  bool hasRowsToSecond = rowsToSecond;
  if (!hasRowsToSecond)
    CMR_CALL( CMRallocStackArray(cmr, &rowsToSecond, matrix->numRows) );
  bool hasColumnsToSecond = columnsToSecond;
  if (!hasColumnsToSecond)
    CMR_CALL( CMRallocStackArray(cmr, &columnsToSecond, matrix->numColumns) );

  char* denseColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseColumn, matrix->numColumns) );

  /* Number of rows of A. */
  size_t numRows = 1;
  size_t extraRow = SIZE_MAX;
  for (size_t row = 0; row < sepa->numRows; ++row)
  {
    if ((sepa->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
    {
      rowsToSecond[row] = numRows;
      secondRowsOrigin[numRows++] = row;
    }
    else
    {
      rowsToSecond[row] = SIZE_MAX;
      if (extraRow == SIZE_MAX && sepa->rowsFlags[row] & CMR_SEPA_MASK_EXTRA)
        extraRow = row;
    }
  }

  /* Number of columns of D. */
  size_t numColumns = 2;
  size_t extraColumn = SIZE_MAX;
  for (size_t column = 0; column < sepa->numColumns; ++column)
  {
    if ((sepa->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
    {
      columnsToSecond[column] = numColumns;
      secondColumnsOrigin[numColumns++] = column;
    }
    else
    {
      columnsToSecond[column] = SIZE_MAX;
      if (extraColumn == SIZE_MAX && sepa->columnsFlags[column] & CMR_SEPA_MASK_EXTRA)
        extraColumn = column;
    }
  }

  assert(extraRow < SIZE_MAX);
  assert(extraColumn < SIZE_MAX);

  secondRowsOrigin[0] = extraRow;
  secondColumnsOrigin[0] = extraColumn;
  secondColumnsOrigin[1] = extraColumn;

  /* Count number of nonzeros, copy column vector and check whether we need to negate a row/column. */
  size_t numNonzeros = 1;
  char scaleTopRight = 0;
  char scaleBottomLeft = 0;
  for (size_t row1 = 0; row1 < numRows; ++row1)
  {
    size_t row = secondRowsOrigin[row1];
    denseColumn[row1] = 0;
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t column1 = columnsToSecond[column];
      if (column1 < SIZE_MAX)
      {
        numNonzeros++;
        if (row == extraRow && scaleTopRight == 0)
        {
          CMRdbgMsg(2, "Scaling entry for top-right is %d, taken from row r%zu, c%zu.\n", matrix->entryValues[e], row+1,
            column+1);
          scaleTopRight = matrix->entryValues[e];
        }
      }
      else if (column == extraColumn && row != extraRow)
      {
        numNonzeros += 2;
        denseColumn[row1] = matrix->entryValues[e];
        if (scaleBottomLeft == 0)
        {
          CMRdbgMsg(2, "Scaling entry for bottom-left is %d, taken from row r%zu, c%zu.\n", matrix->entryValues[e],
            row+1, column+1);
          scaleBottomLeft = matrix->entryValues[e];
        }
      }
    }
  }

  /* Copy the matrix entries. */
  CMR_CALL( CMRchrmatCreate(cmr, psecond, numRows, numColumns, numNonzeros) );
  CMR_CHRMAT* second = *psecond;
  second->numNonzeros = 0;
  for (size_t row1 = 0; row1 < numRows; ++row1)
  {
    second->rowSlice[row1] = second->numNonzeros;
    size_t row = secondRowsOrigin[row1];

    if (row != extraRow)
    {
      /* Nonzero of d. */
      if (denseColumn[row1])
      {
        second->entryColumns[second->numNonzeros] = 0;
        second->entryValues[second->numNonzeros++] = denseColumn[row1] * scaleBottomLeft;
        second->entryColumns[second->numNonzeros] = 1;
        second->entryValues[second->numNonzeros++] = denseColumn[row1] * scaleBottomLeft;
      }
    }
    else
    {
      /* Nonzero above d. */
      second->entryColumns[second->numNonzeros] = 0;
      second->entryValues[second->numNonzeros++] = epsilon;
    }

    /* Row from D and potentially b^T. */
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = matrix->rowSlice[row]; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      size_t column1 = columnsToSecond[column];
      if (column1 < SIZE_MAX)
      {
        second->entryColumns[second->numNonzeros] = column1;
        char scale = (row != extraRow) ? 1 : scaleTopRight;
        second->entryValues[second->numNonzeros++] = matrix->entryValues[e] * scale;
      }
    }
  }
  second->rowSlice[second->numRows] = second->numNonzeros;

  /* Set markers to last row and last two columns. */
  if (secondSpecialRows)
    secondSpecialRows[0] = 0;
  if (secondSpecialColumns)
  {
    secondSpecialColumns[0] = 0;
    secondSpecialColumns[1] = 1;
  }

  /* Free local arrays. */
  CMR_CALL( CMRfreeStackArray(cmr, &denseColumn) );
  if (hasColumnsToSecond)
    CMR_CALL( CMRfreeStackArray(cmr, &columnsToSecond) );
  if (hasRowsToSecond)
    CMR_CALL( CMRfreeStackArray(cmr, &rowsToSecond) );
  if (hasSecondColumnsOrigin)
    CMR_CALL( CMRfreeStackArray(cmr, &secondColumnsOrigin) );
  if (hasSecondRowsOrigin)
    CMR_CALL( CMRfreeStackArray(cmr, &secondRowsOrigin) );

  return CMR_OKAY;
}

CMR_ERROR CMRthreeSumTruemperCompose(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, size_t* firstSpecialRows,
  size_t* firstSpecialColumns, size_t* secondSpecialRows, size_t* secondSpecialColumns, int8_t characteristic,
  CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(first);
  assert(second);
  assert(presult);

  CMRdbgMsg(0, "CMRthreeSumTruemperCompose for a %zux%zu and a %zux%zu matrix.\n", first->numRows, first->numColumns,
    second->numRows, second->numColumns);

  if (!firstSpecialRows || (firstSpecialRows[0] >= first->numRows)
    || (firstSpecialRows[1] >= first->numRows))
  {
    return CMR_ERROR_INPUT;
  }
  if (!firstSpecialColumns || (firstSpecialColumns[0] >= first->numColumns)
    || (firstSpecialColumns[1] >= first->numColumns) || (firstSpecialColumns[2] >= first->numColumns))
  {
    return CMR_ERROR_INPUT;
  }
  if (!secondSpecialRows || (secondSpecialRows[0] >= second->numRows)
    || (secondSpecialRows[1] >= second->numRows) || (secondSpecialRows[2] >= second->numRows))
  {
    return CMR_ERROR_INPUT;
  }
  if (!secondSpecialColumns || (secondSpecialColumns[0] >= second->numColumns)
    || (secondSpecialColumns[1] >= second->numColumns))
  {
    return CMR_ERROR_INPUT;
  }

  CMR_ERROR error = CMR_OKAY;



  return error;
}


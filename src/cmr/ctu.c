#include <cmr/ctu.h>

#include "env_internal.h"

#include <assert.h>
#include <stdint.h>

CMR_ERROR CMRcomplementRowColumn(CMR* cmr, CMR_CHRMAT* matrix, size_t complementRow, size_t complementColumn,
  CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(complementRow < matrix->numRows || complementRow == SIZE_MAX);
  assert(complementColumn < matrix->numColumns || complementColumn == SIZE_MAX);
  assert(presult);

  /* We copy directly if no complementing is requested. */

  if (complementRow == SIZE_MAX && complementColumn == SIZE_MAX)
  {
    CMR_CALL( CMRchrmatCopy(cmr, matrix, presult) );

    return CMR_OKAY;
  }

  /* Create dense copy of matrix because of complementing. */

  bool* dense = NULL;
  size_t size = matrix->numRows * matrix->numColumns;
  CMR_CALL( CMRallocStackArray(cmr, &dense, size) );
  for (size_t i = 0; i < size; ++i)
    dense[i] = false;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowStarts[row];
    size_t beyond = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      if (matrix->entryValues[e] != 1)
        return CMR_ERROR_INPUT;
      dense[matrix->numRows * row + column] = true;
    }
  }

  /* Copy of row to be complemented. */
  char* rows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rows, matrix->numRows) );
  if (complementColumn < SIZE_MAX)
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
      rows[row] = dense[matrix->numRows * row + complementColumn] ? 1 : 0;
  }
  else
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
      rows[row] = 0;
  }

  /* Copy of column to be complemented. */
  char* columns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columns, matrix->numColumns) );
  if (complementRow < SIZE_MAX)
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columns[column] = dense[matrix->numRows * complementRow + column] ? 1 : 0;
  }
  else
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columns[column] = 0;
  }

  size_t countNonzeros = 0;
  char complementRowColumn1 = (complementRow < SIZE_MAX && complementColumn < SIZE_MAX
    && dense[complementRow * matrix->numRows + complementColumn]) ? 1 : 0;
  size_t entry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
    {
      if (row == complementRow)
      {
        if (column != complementColumn && complementRowColumn1)
          dense[entry] = !dense[entry];
      }
      else
      {
        if (column == complementColumn)
        {
          if (complementRowColumn1)
            dense[entry] = !dense[entry];
        }
        else
        {
          if (complementRowColumn1 + rows[row] + columns[column] % 2 == 1)
            dense[entry] = !dense[entry];
        }
      }
      if (dense[entry])
        ++countNonzeros;
      ++entry;
    }
  }

  /* Copy dense to new sparse matrix. */
  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, countNonzeros) );
  CMR_CHRMAT* result = *presult;
  entry = 0;
  countNonzeros = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    result->rowStarts[row] = countNonzeros;
    for (size_t column = 0; column < matrix->numColumns; ++column)
    {
      if (dense[entry])
      {
        result->entryColumns[countNonzeros] = column;
        result->entryValues[countNonzeros] = 1;
      }
      ++entry;
    }
  }

  CMR_CALL( CMRfreeStackArray(cmr, &columns) );
  CMR_CALL( CMRfreeStackArray(cmr, &rows) );
  CMR_CALL( CMRfreeStackArray(cmr, &dense) );

  return CMR_OKAY;
}

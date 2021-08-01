#include "common.h"

#include <assert.h>
#include <stdlib.h>

#include <stdio.h>

CMR_ERROR stringToDoubleMatrix(CMR* cmr, CMR_DBLMAT** matrix, const char* string)
{
  assert(cmr);
  char* end;
  int numRows, numColumns;

  /* Read size of matrix. */
  numRows = strtol(string, &end, 10);
  assert(end > string);
  string = end;
  numColumns = strtol(string, &end, 10);
  assert(end > string);
  string = end;

  CMRdblmatCreate(cmr, matrix, numRows, numColumns, numRows * numColumns);

  (*matrix)->numNonzeros = 0;
  for (int row = 0; row < (*matrix)->numRows; ++row)
  {
    (*matrix)->rowStarts[row] = (*matrix)->numNonzeros;
    for (int column = 0; column < (*matrix)->numColumns; ++column)
    {
      double x = strtod(string, &end);
      assert(end > string);
      string = end;

      if (x != 0.0)
      {
        (*matrix)->entryColumns[(*matrix)->numNonzeros] = column;
        (*matrix)->entryValues[(*matrix)->numNonzeros] = x;
        (*matrix)->numNonzeros++;
      }
    }
  }

  return CMR_OKAY;
}

CMR_ERROR stringToIntMatrix(CMR* cmr, CMR_INTMAT** matrix, const char* string)
{
  assert(cmr);
  char* end;
  int numRows, numColumns;

  /* Read size of matrix. */
  numRows = strtol(string, &end, 10);
  assert(end > string);
  string = end;
  numColumns = strtol(string, &end, 10);
  assert(end > string);
  string = end;

  CMRintmatCreate(cmr, matrix, numRows, numColumns, numRows * numColumns);

  (*matrix)->numNonzeros = 0;
  for (int row = 0; row < (*matrix)->numRows; ++row)
  {
    (*matrix)->rowStarts[row] = (*matrix)->numNonzeros;
    for (int column = 0; column < (*matrix)->numColumns; ++column)
    {
      int x = strtol(string, &end, 10);
      assert(end > string);
      string = end;

      if (x != 0)
      {
        (*matrix)->entryColumns[(*matrix)->numNonzeros] = column;
        (*matrix)->entryValues[(*matrix)->numNonzeros] = x;
        (*matrix)->numNonzeros++;
      }
    }
  }

  return CMR_OKAY;
}

CMR_ERROR stringToCharMatrix(CMR* cmr, CMR_CHRMAT** matrix, const char* string)
{
  assert(cmr);
  char* end;
  int numRows, numColumns;

  /* Read size of matrix. */
  numRows = strtol(string, &end, 10);
  assert(end > string);
  string = end;
  numColumns = strtol(string, &end, 10);
  assert(end > string);
  string = end;

  CMR_CALL( CMRchrmatCreate(cmr, matrix, numRows, numColumns, numRows * numColumns) );

  (*matrix)->numNonzeros = 0;
  for (int row = 0; row < (*matrix)->numRows; ++row)
  {
    (*matrix)->rowStarts[row] = (*matrix)->numNonzeros;
    for (int column = 0; column < (*matrix)->numColumns; ++column)
    {
      int x = strtol(string, &end, 10);
      if (end == string)
        return CMR_ERROR_INPUT;
      string = end;

      if (x != 0)
      {
        (*matrix)->entryColumns[(*matrix)->numNonzeros] = column;
        (*matrix)->entryValues[(*matrix)->numNonzeros] = x;
        (*matrix)->numNonzeros++;
      }
    }
  }

  return CMR_OKAY;
}

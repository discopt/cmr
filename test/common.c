#include "common.h"

#include <assert.h>
#include <stdlib.h>

#include <stdio.h>

CMR_ERROR stringToDoubleMatrix(CMR* cmr, CMR_DBLMAT** pmatrix, const char* string)
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

  CMRdblmatCreate(cmr, pmatrix, numRows, numColumns, numRows * numColumns);

  (*pmatrix)->numNonzeros = 0;
  for (size_t row = 0; row < (*pmatrix)->numRows; ++row)
  {
    (*pmatrix)->rowSlice[row] = (*pmatrix)->numNonzeros;
    for (size_t column = 0; column < (*pmatrix)->numColumns; ++column)
    {
      double x = strtod(string, &end);
      assert(end > string);
      string = end;

      if (x != 0.0)
      {
        (*pmatrix)->entryColumns[(*pmatrix)->numNonzeros] = column;
        (*pmatrix)->entryValues[(*pmatrix)->numNonzeros] = x;
        (*pmatrix)->numNonzeros++;
      }
    }
  }
  (*pmatrix)->rowSlice[numRows] = (*pmatrix)->numNonzeros;

  return CMR_OKAY;
}

CMR_ERROR stringToIntMatrix(CMR* cmr, CMR_INTMAT** pmatrix, const char* string)
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

  CMRintmatCreate(cmr, pmatrix, numRows, numColumns, numRows * numColumns);

  (*pmatrix)->numNonzeros = 0;
  for (size_t row = 0; row < (*pmatrix)->numRows; ++row)
  {
    (*pmatrix)->rowSlice[row] = (*pmatrix)->numNonzeros;
    for (size_t column = 0; column < (*pmatrix)->numColumns; ++column)
    {
      int x = strtol(string, &end, 10);
      assert(end > string);
      string = end;

      if (x != 0)
      {
        (*pmatrix)->entryColumns[(*pmatrix)->numNonzeros] = column;
        (*pmatrix)->entryValues[(*pmatrix)->numNonzeros] = x;
        (*pmatrix)->numNonzeros++;
      }
    }
  }
  (*pmatrix)->rowSlice[numRows] = (*pmatrix)->numNonzeros;

  return CMR_OKAY;
}

CMR_ERROR stringToCharMatrix(CMR* cmr, CMR_CHRMAT** pmatrix, const char* string)
{
  assert(cmr);
  char* end;
  size_t numRows, numColumns;

  /* Read size of matrix. */
  numRows = strtol(string, &end, 10);
  assert(end > string);
  string = end;
  numColumns = strtol(string, &end, 10);
  assert(end > string);
  string = end;

  CMR_CALL( CMRchrmatCreate(cmr, pmatrix, numRows, numColumns, numRows * numColumns) );

  (*pmatrix)->numNonzeros = 0;
  for (size_t row = 0; row < (*pmatrix)->numRows; ++row)
  {
    (*pmatrix)->rowSlice[row] = (*pmatrix)->numNonzeros;
    for (size_t column = 0; column < (*pmatrix)->numColumns; ++column)
    {
      int x = strtol(string, &end, 10);
      if (end == string)
        return CMR_ERROR_INPUT;
      string = end;

      if (x != 0)
      {
        (*pmatrix)->entryColumns[(*pmatrix)->numNonzeros] = column;
        (*pmatrix)->entryValues[(*pmatrix)->numNonzeros] = x;
        (*pmatrix)->numNonzeros++;
      }
    }
  }
  (*pmatrix)->rowSlice[numRows] = (*pmatrix)->numNonzeros;

  return CMR_OKAY;
}

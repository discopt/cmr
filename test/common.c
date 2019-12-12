#include "common.h"

#include <assert.h>
#include <stdlib.h>

#include <stdio.h>

TU_SPARSE_DOUBLE stringToSparseDouble(const char* string)
{
  TU_SPARSE_DOUBLE matrix;
  char* end;
  int maxNonzeros = 256;

  /* Read size of matrix. */
  matrix.numRows = strtol(string, &end, 10);
  assert(end > string);
  string = end;
  matrix.numColumns = strtol(string, &end, 10);
  assert(end > string);
  string = end;

  matrix.majorStarts = (int*) malloc(matrix.numRows * sizeof(int));
  matrix.numNonzeros = 0;
  matrix.entryMinors = (int*) malloc(maxNonzeros * sizeof(int));
  matrix.entryValues = (double*) malloc(maxNonzeros * sizeof(double));

  for (int row = 0; row < matrix.numRows; ++row)
  {
    matrix.majorStarts[row] = matrix.numNonzeros;
    for (int column = 0; column < matrix.numColumns; ++column)
    {
      double x = strtod(string, &end);
      assert(end > string);
      string = end;

      if (x != 0.0)
      {
        /* Enlarge array if necessary. */
        if (matrix.numNonzeros == maxNonzeros)
        {
          maxNonzeros *= 2;
          matrix.entryMinors = (int*) realloc(matrix.entryMinors, maxNonzeros * sizeof(int));
          matrix.entryValues = (double*) realloc(matrix.entryValues, maxNonzeros * sizeof(double));
        }

        matrix.entryMinors[matrix.numNonzeros] = column;
        matrix.entryValues[matrix.numNonzeros] = x;
        matrix.numNonzeros++;
      }
    }
  }

  return matrix;
}

TU_SPARSE_INT stringToSparseInt(const char* string)
{
  TU_SPARSE_INT matrix;
  char* end;
  int maxNonzeros = 256;

  /* Read size of matrix. */
  matrix.numRows = strtol(string, &end, 10);
  assert(end > string);
  string = end;
  matrix.numColumns = strtol(string, &end, 10);
  assert(end > string);
  string = end;

  matrix.rowStarts = (int*) malloc(matrix.numRows * sizeof(int));
  matrix.numNonzeros = 0;
  matrix.entryColumns = (int*) malloc(maxNonzeros * sizeof(int));
  matrix.entryValues = (int*) malloc(maxNonzeros * sizeof(int));

  for (int row = 0; row < matrix.numRows; ++row)
  {
    matrix.rowStarts[row] = matrix.numNonzeros;
    for (int column = 0; column < matrix.numColumns; ++column)
    {
      int x = strtol(string, &end, 10);
      assert(end > string);
      string = end;

      if (x != 0.0)
      {
        /* Enlarge array if necessary. */
        if (matrix.numNonzeros == maxNonzeros)
        {
          maxNonzeros *= 2;
          matrix.entryColumns = (int*) realloc(matrix.entryColumns, maxNonzeros * sizeof(int));
          matrix.entryValues = (int*) realloc(matrix.entryValues, maxNonzeros * sizeof(int));
        }

        matrix.entryColumns[matrix.numNonzeros] = column;
        matrix.entryValues[matrix.numNonzeros] = x;
        matrix.numNonzeros++;
      }
    }
  }

  return matrix;
}

TU_SPARSE_CHAR stringToSparseChar(const char* string)
{
  TU_SPARSE_CHAR matrix;
  char* end;
  int maxNonzeros = 256;

  /* Read size of matrix. */
  matrix.numRows = strtol(string, &end, 10);
  assert(end > string);
  string = end;
  matrix.numColumns = strtol(string, &end, 10);
  assert(end > string);
  string = end;

  matrix.majorStarts = (int*) malloc(matrix.numRows * sizeof(int));
  matrix.numNonzeros = 0;
  matrix.entryMinors = (int*) malloc(maxNonzeros * sizeof(int));
  matrix.entryValues = (char*) malloc(maxNonzeros * sizeof(char));

  for (int row = 0; row < matrix.numRows; ++row)
  {
    matrix.majorStarts[row] = matrix.numNonzeros;
    for (int column = 0; column < matrix.numColumns; ++column)
    {
      char x = strtol(string, &end, 10);
      assert(end > string);
      string = end;

      if (x != 0.0)
      {
        /* Enlarge array if necessary. */
        if (matrix.numNonzeros == maxNonzeros)
        {
          maxNonzeros *= 2;
          matrix.entryMinors = (int*) realloc(matrix.entryMinors, maxNonzeros * sizeof(int));
          matrix.entryValues = (char*) realloc(matrix.entryValues, maxNonzeros * sizeof(char));
        }

        matrix.entryMinors[matrix.numNonzeros] = column;
        matrix.entryValues[matrix.numNonzeros] = x;
        matrix.numNonzeros++;
      }
    }
  }

  return matrix;
}

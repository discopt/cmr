#include <tu/matrix.h>

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "sort.h"
#include "env_internal.h"

TU_ERROR TUdblmatCreate(TU* tu, TU_DBLMAT** matrix, int numRows, int numColumns,
  int numNonzeros)
{
  assert(matrix);
  assert(*matrix == NULL);

  TU_CALL( TUallocBlock(tu, matrix) );
  (*matrix)->numRows = numRows;
  (*matrix)->numColumns = numColumns;
  (*matrix)->numNonzeros = numNonzeros;
  (*matrix)->rowStarts = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  TU_CALL( TUallocBlockArray(tu, &(*matrix)->rowStarts, numRows + 1) );
  if (numNonzeros > 0)
  {
    TU_CALL( TUallocBlockArray(tu, &(*matrix)->entryColumns, numNonzeros) );
    TU_CALL( TUallocBlockArray(tu, &(*matrix)->entryValues, numNonzeros) );
  }

  return TU_OKAY;
}

TU_ERROR TUdblmatFree(TU* tu, TU_DBLMAT** matrix)
{
  assert(matrix);
  assert(*matrix);
  assert((*matrix)->rowStarts);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryColumns);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryValues);

  TU_CALL( TUfreeBlockArray(tu, &(*matrix)->rowStarts) );
  if ((*matrix)->numNonzeros > 0)
  {
    TU_CALL( TUfreeBlockArray(tu, &(*matrix)->entryColumns) );
    TU_CALL( TUfreeBlockArray(tu, &(*matrix)->entryValues) );
  }
  TU_CALL( TUfreeBlock(tu, matrix) );

  return TU_OKAY;
}

TU_ERROR TUdblmatChangeNumNonzeros(TU* tu, TU_DBLMAT* matrix, int newNumNonzeros)
{
  assert(tu);
  assert(matrix);

  TU_CALL( TUreallocBlockArray(tu, &matrix->entryColumns, newNumNonzeros) );
  TU_CALL( TUreallocBlockArray(tu, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return TU_OKAY;
}

TU_ERROR TUdblmatCopy(TU* tu, TU_DBLMAT* matrix, TU_DBLMAT** result)
{
  assert(tu);
  assert(matrix);
  assert(result);
  assert(*result == NULL);

  TU_CALL( TUdblmatCreate(tu, result, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  for (int row = 0; row < matrix->numRows; ++row)
    (*result)->rowStarts[row] = matrix->rowStarts[row];
  (*result)->rowStarts[matrix->numRows] = matrix->numNonzeros;
  for (int entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    (*result)->entryColumns[entry] = matrix->entryColumns[entry];
    (*result)->entryValues[entry] = matrix->entryValues[entry];
  }

  return TU_OKAY;
}

TU_ERROR TUdblmatTranspose(TU* tu, TU_DBLMAT* matrix, TU_DBLMAT** result)
{
  assert(tu);
  assert(matrix);
  assert(result);
  assert(*result == NULL);
  assert(TUdblmatCheckSorted(matrix));

  TU_CALL( TUdblmatCreate(tu, result, matrix->numColumns, matrix->numRows, matrix->numNonzeros) );

  /* Count number of nonzeros in each column, storing in the next entry. */
  for (int c = 0; c <= matrix->numColumns; ++c)
    (*result)->rowStarts[c] = 0;
  for (int e = 0; e < matrix->numNonzeros; ++e)
    (*result)->rowStarts[matrix->entryColumns[e] + 1]++;

  /* Compute start indices for columns. */
  for (int c = 1; c < matrix->numColumns; ++c)
    (*result)->rowStarts[c] += (*result)->rowStarts[c-1];

  /* Create nonzeros. */
  for (int row = 0; row < matrix->numRows; ++row)
  {
    int begin = matrix->rowStarts[row];
    int end = row + 1 < matrix->numRows ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      int column = matrix->entryColumns[entry];
      int transEntry = (*result)->rowStarts[column];
      (*result)->entryColumns[transEntry] = row;
      (*result)->entryValues[transEntry] = matrix->entryValues[entry];
      (*result)->rowStarts[column]++;
    }
  }

  /* We shifted rowStarts of *result, so we shift it back. */
  for (int c = matrix->numColumns; c > 0; --c)
    (*result)->rowStarts[c] = (*result)->rowStarts[c-1];
  (*result)->rowStarts[0] = 0;

  return TU_OKAY;
}

TU_ERROR TUintmatCreate(TU* tu, TU_INTMAT** matrix, int numRows, int numColumns, int numNonzeros)
{
  assert(matrix);
  assert(*matrix == NULL);

  TU_CALL( TUallocBlock(tu, matrix) );
  (*matrix)->numRows = numRows;
  (*matrix)->numColumns = numColumns;
  (*matrix)->numNonzeros = numNonzeros;
  (*matrix)->rowStarts = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  TU_CALL( TUallocBlockArray(tu, &(*matrix)->rowStarts, numRows + 1) );
  if (numNonzeros > 0)
  {
    TU_CALL( TUallocBlockArray(tu, &(*matrix)->entryColumns, numNonzeros) );
    TU_CALL( TUallocBlockArray(tu, &(*matrix)->entryValues, numNonzeros) );
  }

  return TU_OKAY;
}

TU_ERROR TUintmatFree(TU* tu, TU_INTMAT** matrix)
{
  assert(matrix);
  assert(*matrix);
  assert((*matrix)->rowStarts);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryColumns);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryValues);

  TU_CALL( TUfreeBlockArray(tu, &(*matrix)->rowStarts) );
  if ((*matrix)->numNonzeros > 0)
  {
    TU_CALL( TUfreeBlockArray(tu, &(*matrix)->entryColumns) );
    TU_CALL( TUfreeBlockArray(tu, &(*matrix)->entryValues) );
  }
  TU_CALL( TUfreeBlock(tu, matrix) );

  return TU_OKAY;
}

TU_ERROR TUintmatChangeNumNonzeros(TU* tu, TU_INTMAT* matrix, int newNumNonzeros)
{
  assert(tu);
  assert(matrix);

  TU_CALL( TUreallocBlockArray(tu, &matrix->entryColumns, newNumNonzeros) );
  TU_CALL( TUreallocBlockArray(tu, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return TU_OKAY;
}

TU_ERROR TUintmatCopy(TU* tu, TU_INTMAT* matrix, TU_INTMAT** result)
{
  assert(tu);
  assert(matrix);
  assert(result);
  assert(*result == NULL);

  TU_CALL( TUintmatCreate(tu, result, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  for (int row = 0; row < matrix->numRows; ++row)
    (*result)->rowStarts[row] = matrix->rowStarts[row];
  (*result)->rowStarts[matrix->numRows] = matrix->numNonzeros;
  for (int entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    (*result)->entryColumns[entry] = matrix->entryColumns[entry];
    (*result)->entryValues[entry] = matrix->entryValues[entry];
  }

  return TU_OKAY;
}

TU_ERROR TUchrmatCreate(TU* tu, TU_CHRMAT** matrix, int numRows, int numColumns, int numNonzeros)
{
  assert(matrix);
  assert(*matrix == NULL);

  TU_CALL( TUallocBlock(tu, matrix) );
  (*matrix)->numRows = numRows;
  (*matrix)->numColumns = numColumns;
  (*matrix)->numNonzeros = numNonzeros;
  (*matrix)->rowStarts = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  TU_CALL( TUallocBlockArray(tu, &(*matrix)->rowStarts, numRows + 1) );
  if (numNonzeros > 0)
  {
    TU_CALL( TUallocBlockArray(tu, &(*matrix)->entryColumns, numNonzeros) );
    TU_CALL( TUallocBlockArray(tu, &(*matrix)->entryValues, numNonzeros) );
  }

  return TU_OKAY;
}

TU_ERROR TUchrmatFree(TU* tu, TU_CHRMAT** matrix)
{
  assert(matrix);
  assert(*matrix);
  assert((*matrix)->rowStarts);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryColumns);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryValues);

  TU_CALL( TUfreeBlockArray(tu, &(*matrix)->rowStarts) );
  if ((*matrix)->numNonzeros > 0)
  {
    TU_CALL( TUfreeBlockArray(tu, &(*matrix)->entryColumns) );
    TU_CALL( TUfreeBlockArray(tu, &(*matrix)->entryValues) );
  }
  TU_CALL( TUfreeBlock(tu, matrix) );

  return TU_OKAY;
}

TU_ERROR TUchrmatChangeNumNonzeros(TU* tu, TU_CHRMAT* matrix, int newNumNonzeros)
{
  assert(tu);
  assert(matrix);

  TU_CALL( TUreallocBlockArray(tu, &matrix->entryColumns, newNumNonzeros) );
  TU_CALL( TUreallocBlockArray(tu, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return TU_OKAY;
}

TU_ERROR TUchrmatCopy(TU* tu, TU_CHRMAT* matrix, TU_CHRMAT** result)
{
  assert(tu);
  assert(matrix);
  assert(result);
  assert(*result == NULL);

  TU_CALL( TUchrmatCreate(tu, result, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  for (int row = 0; row < matrix->numRows; ++row)
    (*result)->rowStarts[row] = matrix->rowStarts[row];
  (*result)->rowStarts[matrix->numRows] = matrix->numNonzeros;
  for (int entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    (*result)->entryColumns[entry] = matrix->entryColumns[entry];
    (*result)->entryValues[entry] = matrix->entryValues[entry];
  }

  return TU_OKAY;
}


TU_ERROR TUchrmatTranspose(TU* tu, TU_CHRMAT* matrix, TU_CHRMAT** result)
{
  assert(tu);
  assert(matrix);
  assert(result);
  assert(*result == NULL);
  assert(TUchrmatCheckSorted(matrix));

  TU_CALL( TUchrmatCreate(tu, result, matrix->numColumns, matrix->numRows, matrix->numNonzeros) );

  /* Count number of nonzeros in each column, storing in the next entry. */
  for (int c = 0; c <= matrix->numColumns; ++c)
    (*result)->rowStarts[c] = 0;
  for (int e = 0; e < matrix->numNonzeros; ++e)
    (*result)->rowStarts[matrix->entryColumns[e] + 1]++;

  /* Compute start indices for columns. */
  for (int c = 1; c < matrix->numColumns; ++c)
    (*result)->rowStarts[c] += (*result)->rowStarts[c-1];

  /* Create nonzeros. */
  for (int row = 0; row < matrix->numRows; ++row)
  {
    int begin = matrix->rowStarts[row];
    int end = row + 1 < matrix->numRows ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      int column = matrix->entryColumns[entry];
      int transEntry = (*result)->rowStarts[column];
      (*result)->entryColumns[transEntry] = row;
      (*result)->entryValues[transEntry] = matrix->entryValues[entry];
      (*result)->rowStarts[column]++;
    }
  }

  /* We shifted rowStarts of *result, so we shift it back. */
  for (int c = matrix->numColumns; c > 0; --c)
    (*result)->rowStarts[c] = (*result)->rowStarts[c-1];
  (*result)->rowStarts[0] = 0;

  return TU_OKAY;
}

TU_ERROR TUdblmatPrintSparse(FILE* stream, TU_DBLMAT* matrix)
{
  assert(stream);
  assert(matrix);

  fprintf(stream, "%d %d %d\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (int row = 0; row < matrix->numRows; ++row)
  {
    int first = matrix->rowStarts[row];
    int beyond = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int entry = first; entry < beyond; ++entry)
      fprintf(stream, "%d %d  %f\n", row, matrix->entryColumns[entry], matrix->entryValues[entry]);
  }

  return TU_OKAY;
}

TU_ERROR TUintmatPrintSparse(FILE* stream, TU_INTMAT* matrix)
{
  assert(stream);
  assert(matrix);

  fprintf(stream, "%d %d %d\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (int row = 0; row < matrix->numRows; ++row)
  {
    int first = matrix->rowStarts[row];
    int beyond = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int entry = first; entry < beyond; ++entry)
      fprintf(stream, "%d %d  %d\n", row, matrix->entryColumns[entry], matrix->entryValues[entry]);
  }

  return TU_OKAY;
}

TU_ERROR TUchrmatPrintSparse(FILE* stream, TU_CHRMAT* matrix)
{
  assert(stream);
  assert(matrix);

  fprintf(stream, "%d %d %d\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (int row = 0; row < matrix->numRows; ++row)
  {
    int first = matrix->rowStarts[row];
    int beyond = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int entry = first; entry < beyond; ++entry)
      fprintf(stream, "%d %d  %d\n", row, matrix->entryColumns[entry], matrix->entryValues[entry]);
  }

  return TU_OKAY;
}



TU_ERROR TUdblmatPrintDense(FILE* stream, TU_DBLMAT* matrix, char zeroChar, bool header)
{
  assert(stream != NULL);
  assert(matrix != NULL);
  double* rowEntries = (double*) calloc(matrix->numColumns, sizeof(double));

  fprintf(stream, "%d %d\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (int column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%d ", column % 10);
    fputs("\n  ", stream);
    for (int column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (int row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%d| ", row % 10);
    int start = matrix->rowStarts[row];
    int end = row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;
    for (int i = start; i < end; ++i)
      rowEntries[matrix->entryColumns[i]] = matrix->entryValues[i];
    for (int column = 0; column < matrix->numColumns; ++column)
    {
      double x = rowEntries[column];
      if (x == 0.0)
        fprintf(stream, "%c ", zeroChar);
      else
        fprintf(stream, "%f ", x);
    }
    for (int i = start; i < end; ++i)
      rowEntries[matrix->entryColumns[i]] = 0.0;
    fputc('\n', stream);
  }

  free(rowEntries);

  return TU_OKAY;
}

TU_ERROR TUintmatPrintDense(FILE* stream, TU_INTMAT* matrix, char zeroChar, bool header)
{
  assert(stream != NULL);
  assert(matrix != NULL);
  int* rowEntries = (int*) calloc(matrix->numColumns, sizeof(int));

  fprintf(stream, "%d %d\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (int column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%d ", column % 10);
    fputs("\n  ", stream);
    for (int column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (int row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%d| ", row % 10);
    int start = matrix->rowStarts[row];
    int end = row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;
    for (int i = start; i < end; ++i)
      rowEntries[matrix->entryColumns[i]] = matrix->entryValues[i];
    for (int column = 0; column < matrix->numColumns; ++column)
    {
      int x = rowEntries[column];
      if (x == 0.0)
        fprintf(stream, "%c ", zeroChar);
      else
        fprintf(stream, "%d ", x);
    }
    for (int i = start; i < end; ++i)
      rowEntries[matrix->entryColumns[i]] = 0.0;
    fputc('\n', stream);
  }

  free(rowEntries);

  return TU_OKAY;
}

TU_ERROR TUchrmatPrintDense(FILE* stream, TU_CHRMAT* matrix, char zeroChar, bool header)
{
  assert(stream != NULL);
  assert(matrix != NULL);
  char* rowEntries = (char*) calloc(matrix->numColumns, sizeof(char));

  fprintf(stream, "%d %d\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (int column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%d ", column % 10);
    fputs("\n  ", stream);
    for (int column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (int row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%d| ", row % 10);
    int start = matrix->rowStarts[row];
    int end = row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;
    for (int i = start; i < end; ++i)
      rowEntries[matrix->entryColumns[i]] = matrix->entryValues[i];
    for (int column = 0; column < matrix->numColumns; ++column)
    {
      char x = rowEntries[column];
      if (x == 0.0)
        fprintf(stream, "%c ", zeroChar);
      else
        fprintf(stream, "%d ", x);
    }
    for (int i = start; i < end; ++i)
      rowEntries[matrix->entryColumns[i]] = 0.0;
    fputc('\n', stream);
  }

  free(rowEntries);

  return TU_OKAY;
}

typedef struct
{
  int row;
  int column;
  double value;
} DblNonzero;

int compareDblNonzeros(const void* pa, const void* pb)
{
  int aRow = ((const DblNonzero*)pa)->row;
  int bRow = ((const DblNonzero*)pb)->row;
  if (aRow < bRow)
    return -1;
  else if (aRow > bRow)
    return +1;

  int aColumn = ((const DblNonzero*)pa)->column;
  int bColumn = ((const DblNonzero*)pb)->column;
  return aColumn - bColumn;
}

TU_ERROR TUdblmatCreateFromSparseStream(TU* tu, TU_DBLMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%d %d %d", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return TU_ERROR_INPUT;

  /* Read all nonzeros. */

  DblNonzero* nonzeros = NULL;
  TU_CALL( TUallocStackArray(tu, &nonzeros, numNonzeros) );
  int entry = 0;
  for (int i = 0; i < numNonzeros; ++i)
  {
    int row;
    int column;
    double value;
    numRead = fscanf(stream, "%d %d %lf", &row, &column, &value);
    if (numRead < 3 || row < 0 || column < 0 || row >= numRows || column >= numColumns)
    {
      TU_CALL( TUfreeStackArray(tu, &nonzeros) );
      return TU_ERROR_INPUT;
    }
    if (value != 0.0)
    {
      nonzeros[entry].row = row;
      nonzeros[entry].column = column;
      nonzeros[entry].value = value;
      ++entry;
    }
  }
  numNonzeros = entry;

  /* We sort all nonzeros by row and then by column. */
  TU_CALL( TUsort(tu, numNonzeros, nonzeros, sizeof(DblNonzero), compareDblNonzeros) );

  TU_CALL( TUdblmatCreate(tu, pmatrix, numRows, numColumns, numNonzeros) );
  int previousRow = -1;
  int previousColumn = -1;
  int* pentryColumn = (*pmatrix)->entryColumns;
  double* pentryValue = (*pmatrix)->entryValues;
  for (int entry = 0; entry < numNonzeros; ++entry)
  {
    int row = nonzeros[entry].row;
    int column = nonzeros[entry].column;
    if (row == previousRow && column == previousColumn)
    {
      TU_CALL( TUfreeStackArray(tu, &nonzeros) );
      TU_CALL( TUdblmatFree(tu, pmatrix) );
      return TU_ERROR_INPUT;
    }
    while (previousRow < row)
    {
      ++previousRow;
      (*pmatrix)->rowStarts[previousRow] = entry;
    }
    *pentryColumn++ = nonzeros[entry].column;
    *pentryValue++ = nonzeros[entry].value;
    previousColumn = column;
  }
  while (previousRow < numRows)
  {
    ++previousRow;
    (*pmatrix)->rowStarts[previousRow] = numNonzeros;
  }

  TU_CALL( TUfreeStackArray(tu, &nonzeros) );

  return TU_OKAY;
}

typedef struct
{
  int row;
  int column;
  int value;
} IntNonzero;

int compareIntNonzeros(const void* pa, const void* pb)
{
  int aRow = ((const IntNonzero*)pa)->row;
  int bRow = ((const IntNonzero*)pb)->row;
  if (aRow < bRow)
    return -1;
  else if (aRow > bRow)
    return +1;

  int aColumn = ((const IntNonzero*)pa)->column;
  int bColumn = ((const IntNonzero*)pb)->column;
  return aColumn - bColumn;
}

TU_ERROR TUintmatCreateFromSparseStream(TU* tu, TU_INTMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%d %d %d", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return TU_ERROR_INPUT;

  /* Read all nonzeros. */

  IntNonzero* nonzeros = NULL;
  TU_CALL( TUallocStackArray(tu, &nonzeros, numNonzeros) );
  int entry = 0;
  for (int i = 0; i < numNonzeros; ++i)
  {
    int row;
    int column;
    int value;
    numRead = fscanf(stream, "%d %d %d", &row, &column, &value);
    if (numRead < 3 || row < 0 || column < 0 || row >= numRows || column >= numColumns)
    {
      TU_CALL( TUfreeStackArray(tu, &nonzeros) );
      return TU_ERROR_INPUT;
    }
    if (value != 0)
    {
      nonzeros[entry].row = row;
      nonzeros[entry].column = column;
      nonzeros[entry].value = value;
      ++entry;
    }
  }
  numNonzeros = entry;

  /* We sort all nonzeros by row and then by column. */
  TU_CALL( TUsort(tu, numNonzeros, nonzeros, sizeof(IntNonzero), compareIntNonzeros) );

  TU_CALL( TUintmatCreate(tu, pmatrix, numRows, numColumns, numNonzeros) );
  int previousRow = -1;
  int previousColumn = -1;
  int* pentryColumn = (*pmatrix)->entryColumns;
  int* pentryValue = (*pmatrix)->entryValues;
  for (int entry = 0; entry < numNonzeros; ++entry)
  {
    int row = nonzeros[entry].row;
    int column = nonzeros[entry].column;
    if (row == previousRow && column == previousColumn)
    {
      TU_CALL( TUfreeStackArray(tu, &nonzeros) );
      TU_CALL( TUintmatFree(tu, pmatrix) );
      return TU_ERROR_INPUT;
    }
    while (previousRow < row)
    {
      ++previousRow;
      (*pmatrix)->rowStarts[previousRow] = entry;
    }
    *pentryColumn++ = nonzeros[entry].column;
    *pentryValue++ = nonzeros[entry].value;
    previousColumn = column;
  }
  while (previousRow < numRows)
  {
    ++previousRow;
    (*pmatrix)->rowStarts[previousRow] = numNonzeros;
  }

  TU_CALL( TUfreeStackArray(tu, &nonzeros) );

  return TU_OKAY;
}

typedef struct
{
  int row;
  int column;
  char value;
} ChrNonzero;

int compareChrNonzeros(const void* pa, const void* pb)
{
  int aRow = ((const ChrNonzero*)pa)->row;
  int bRow = ((const ChrNonzero*)pb)->row;
  if (aRow < bRow)
    return -1;
  else if (aRow > bRow)
    return +1;

  int aColumn = ((const ChrNonzero*)pa)->column;
  int bColumn = ((const ChrNonzero*)pb)->column;
  return aColumn - bColumn;
}

TU_ERROR TUchrmatCreateFromSparseStream(TU* tu, TU_CHRMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%d %d %d", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return TU_ERROR_INPUT;

  /* Read all nonzeros. */

  ChrNonzero* nonzeros = NULL;
  TU_CALL( TUallocStackArray(tu, &nonzeros, numNonzeros) );
  int entry = 0;
  for (int i = 0; i < numNonzeros; ++i)
  {
    int row;
    int column;
    int value;
    numRead = fscanf(stream, "%d %d %d", &row, &column, &value);
    if (numRead < 3 || row < 0 || column < 0 || row >= numRows || column >= numColumns || value < CHAR_MIN
      || value > CHAR_MAX)
    {
      TU_CALL( TUfreeStackArray(tu, &nonzeros) );
      return TU_ERROR_INPUT;
    }
    if (value != 0)
    {
      nonzeros[entry].row = row;
      nonzeros[entry].column = column;
      nonzeros[entry].value = value;
      ++entry;
    }
  }
  numNonzeros = entry;

  /* We sort all nonzeros by row and then by column. */
  TU_CALL( TUsort(tu, numNonzeros, nonzeros, sizeof(ChrNonzero), compareChrNonzeros) );

  TU_CALL( TUchrmatCreate(tu, pmatrix, numRows, numColumns, numNonzeros) );
  int previousRow = -1;
  int previousColumn = -1;
  int* pentryColumn = (*pmatrix)->entryColumns;
  char* pentryValue = (*pmatrix)->entryValues;
  for (int entry = 0; entry < numNonzeros; ++entry)
  {
    int row = nonzeros[entry].row;
    int column = nonzeros[entry].column;
    if (row == previousRow && column == previousColumn)
    {
      TU_CALL( TUfreeStackArray(tu, &nonzeros) );
      TU_CALL( TUchrmatFree(tu, pmatrix) );
      return TU_ERROR_INPUT;
    }
    while (previousRow < row)
    {
      ++previousRow;
      (*pmatrix)->rowStarts[previousRow] = entry;
    }
    *pentryColumn++ = nonzeros[entry].column;
    *pentryValue++ = nonzeros[entry].value;
    previousColumn = column;
  }
  while (previousRow < numRows)
  {
    ++previousRow;
    (*pmatrix)->rowStarts[previousRow] = numNonzeros;
  }

  TU_CALL( TUfreeStackArray(tu, &nonzeros) );

  return TU_OKAY;
}

TU_ERROR TUdblmatCreateFromDenseStream(TU* tu, TU_DBLMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns;
  int numRead = fscanf(stream, "%d %d", &numRows, &numColumns);
  if (numRead < 2)
    return TU_ERROR_INPUT;

  TU_CALL( TUdblmatCreate(tu, pmatrix, numRows, numColumns, 0) );

  /* Initial memory. */
  int memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  int* entryColumns = NULL;
  double* entryValues = NULL;
  TU_CALL( TUallocBlockArray(tu, &entryColumns, memEntries) );
  TU_CALL( TUallocBlockArray(tu, &entryValues, memEntries) );

  int entry = 0;
  for (int row = 0; row < numRows; ++row)
  {
    (*pmatrix)->rowStarts[row] = entry;
    for (int column = 0; column < numColumns; ++column)
    {
      double x;
      numRead = fscanf(stream, "%lf", &x);
      if (numRead < 1)
        return TU_ERROR_INPUT;

      if (x == 0.0)
        continue;

      if (entry == memEntries)
      {
        memEntries = 2*memEntries;
        TU_CALL( TUreallocBlockArray(tu, &entryColumns, memEntries) );
        TU_CALL( TUreallocBlockArray(tu, &entryValues, memEntries) );
      }

      entryColumns[entry] = column;
      entryValues[entry] = x;
      ++entry;
    }
  }
  (*pmatrix)->rowStarts[numRows] = entry;

  /* Make arrays smaller again. */
  if (entry < memEntries)
  {
    TU_CALL( TUreallocBlockArray(tu, &entryColumns, entry) );
    TU_CALL( TUreallocBlockArray(tu, &entryValues, entry) );
  }

  (*pmatrix)->entryColumns = entryColumns;
  (*pmatrix)->entryValues = entryValues;
  (*pmatrix)->numNonzeros = entry;

  return TU_OKAY;
}

TU_ERROR TUintmatCreateFromDenseStream(TU* tu, TU_INTMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns;
  int numRead = fscanf(stream, "%d %d", &numRows, &numColumns);
  if (numRead < 2)
    return TU_ERROR_INPUT;

  TU_CALL( TUintmatCreate(tu, pmatrix, numRows, numColumns, 0) );

  /* Initial memory. */
  int memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  int* entryColumns = NULL;
  int* entryValues = NULL;
  TU_CALL( TUallocBlockArray(tu, &entryColumns, memEntries) );
  TU_CALL( TUallocBlockArray(tu, &entryValues, memEntries) );

  int entry = 0;
  for (int row = 0; row < numRows; ++row)
  {
    (*pmatrix)->rowStarts[row] = entry;
    for (int column = 0; column < numColumns; ++column)
    {
      int x;
      numRead = fscanf(stream, "%d", &x);
      if (numRead < 1)
        return TU_ERROR_INPUT;

      if (x == 0.0)
        continue;

      if (entry == memEntries)
      {
        memEntries = 2*memEntries;
        TU_CALL( TUreallocBlockArray(tu, &entryColumns, memEntries) );
        TU_CALL( TUreallocBlockArray(tu, &entryValues, memEntries) );
      }

      entryColumns[entry] = column;
      entryValues[entry] = x;
      ++entry;
    }
  }
  (*pmatrix)->rowStarts[numRows] = entry;

  /* Make arrays smaller again. */
  if (entry < memEntries)
  {
    TU_CALL( TUreallocBlockArray(tu, &entryColumns, entry) );
    TU_CALL( TUreallocBlockArray(tu, &entryValues, entry) );
  }

  (*pmatrix)->entryColumns = entryColumns;
  (*pmatrix)->entryValues = entryValues;
  (*pmatrix)->numNonzeros = entry;

  return TU_OKAY;
}

TU_ERROR TUchrmatCreateFromDenseStream(TU* tu, TU_CHRMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows = 0;
  int numColumns = 0;
  int numRead = fscanf(stream, "%d %d", &numRows, &numColumns);
  if (numRead < 2)
    return TU_ERROR_INPUT;

  TU_CALL( TUchrmatCreate(tu, pmatrix, numRows, numColumns, 0) );

  /* Initial memory. */
  int memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  int* entryColumns = NULL;
  char* entryValues = NULL;
  TU_CALL( TUallocBlockArray(tu, &entryColumns, memEntries) );
  TU_CALL( TUallocBlockArray(tu, &entryValues, memEntries) );

  int entry = 0;
  for (int row = 0; row < numRows; ++row)
  {
    (*pmatrix)->rowStarts[row] = entry;
    for (int column = 0; column < numColumns; ++column)
    {
      int x;
      numRead = fscanf(stream, "%d", &x);
      if (numRead < 1)
        return TU_ERROR_INPUT;

      if (x == 0.0)
        continue;

      if (entry == memEntries)
      {
        memEntries = 2*memEntries;
        TU_CALL( TUreallocBlockArray(tu, &entryColumns, memEntries) );
        TU_CALL( TUreallocBlockArray(tu, &entryValues, memEntries) );
      }

      entryColumns[entry] = column;
      entryValues[entry] = x;
      ++entry;
    }
  }
  (*pmatrix)->rowStarts[numRows] = entry;

  /* Make arrays smaller again. */
  if (entry < memEntries)
  {
    TU_CALL( TUreallocBlockArray(tu, &entryColumns, entry) );
    TU_CALL( TUreallocBlockArray(tu, &entryValues, entry) );
  }

  (*pmatrix)->entryColumns = entryColumns;
  (*pmatrix)->entryValues = entryValues;
  (*pmatrix)->numNonzeros = entry;

  return TU_OKAY;
}

bool TUdblmatCheckEqual(TU_DBLMAT* matrix1, TU_DBLMAT* matrix2)
{
  assert(TUdblmatCheckSorted(matrix1));
  assert(TUdblmatCheckSorted(matrix2));

  if (matrix1->numRows != matrix2->numRows)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;

  for (int row = 0; row < matrix1->numRows; ++row)
  {
    int start1 = matrix1->rowStarts[row];
    int start2 = matrix2->rowStarts[row];
    if (start1 != start2)
      return false;
    int end1 = row + 1 < matrix1->numRows ? matrix1->rowStarts[row] : matrix1->numNonzeros;
    int end2 = row + 1 < matrix2->numRows ? matrix2->rowStarts[row] : matrix2->numNonzeros;
    if (end1 != end2)
      return false;

    for (int i = start1; i < end1; ++i)
    {
      if (matrix1->entryColumns[i] != matrix2->entryColumns[i])
        return false;
      if (matrix1->entryValues[i] != matrix2->entryValues[i])
        return false;
    }
  }

  return true;
}

bool TUintmatCheckEqual(TU_INTMAT* matrix1, TU_INTMAT* matrix2)
{
  assert(TUintmatCheckSorted(matrix1));
  assert(TUintmatCheckSorted(matrix2));

  if (matrix1->numRows != matrix2->numRows)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;

  for (int row = 0; row < matrix1->numRows; ++row)
  {
    int start1 = matrix1->rowStarts[row];
    int start2 = matrix2->rowStarts[row];
    if (start1 != start2)
      return false;
    int end1 = row + 1 < matrix1->numRows ? matrix1->rowStarts[row] : matrix1->numNonzeros;
    int end2 = row + 1 < matrix2->numRows ? matrix2->rowStarts[row] : matrix2->numNonzeros;
    if (end1 != end2)
      return false;

    for (int i = start1; i < end1; ++i)
    {
      if (matrix1->entryColumns[i] != matrix2->entryColumns[i])
        return false;
      if (matrix1->entryValues[i] != matrix2->entryValues[i])
        return false;
    }
  }

  return true;
}

bool TUchrmatCheckEqual(TU_CHRMAT* matrix1, TU_CHRMAT* matrix2)
{
  assert(TUchrmatCheckSorted(matrix1));
  assert(TUchrmatCheckSorted(matrix2));

  if (matrix1->numRows != matrix2->numRows)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;

  for (int row = 0; row < matrix1->numRows; ++row)
  {
    int start1 = matrix1->rowStarts[row];
    int start2 = matrix2->rowStarts[row];
    if (start1 != start2)
      return false;
    int end1 = row + 1 < matrix1->numRows ? matrix1->rowStarts[row + 1] : matrix1->numNonzeros;
    int end2 = row + 1 < matrix2->numRows ? matrix2->rowStarts[row + 1] : matrix2->numNonzeros;
    if (end1 != end2)
      return false;

    for (int i = start1; i < end1; ++i)
    {
      if (matrix1->entryColumns[i] != matrix2->entryColumns[i])
        return false;
      if (matrix1->entryValues[i] != matrix2->entryValues[i])
        return false;
    }
  }

  return true;
}

bool TUdblmatCheckTranspose(TU_DBLMAT* matrix1, TU_DBLMAT* matrix2)
{
  bool result = true;

  assert(matrix1 != NULL);
  assert(matrix2 != NULL);

  if (matrix1->numRows != matrix2->numColumns)
    return false;
  if (matrix1->numColumns != matrix2->numRows)
    return false;
  if (matrix1->numNonzeros != matrix2->numNonzeros)
    return false;

  int* currentColumnEntries = (int*) malloc(matrix1->numColumns * sizeof(int) );
  for (int column = 0; column < matrix2->numRows; ++column)
    currentColumnEntries[column] = matrix2->rowStarts[column];

  for (int row = 0; row < matrix1->numRows; ++row)
  {
    int begin = matrix1->rowStarts[row];
    int end = row + 1 < matrix1->numRows ? matrix1->rowStarts[row + 1] : matrix1->numNonzeros;
    for (int entry1 = begin; entry1 < end; ++entry1)
    {
      int column = matrix1->entryColumns[entry1];
      int entry2 = currentColumnEntries[column];
      if (matrix2->entryColumns[entry2] != row
        || matrix2->entryValues[entry2] != matrix1->entryValues[entry1])
      {
        result = false;
        goto cleanup;
      }
      currentColumnEntries[column]++;
    }
  }

cleanup:

  free(currentColumnEntries);

  return result;
}

bool TUintmatCheckTranspose(TU_INTMAT* matrix1, TU_INTMAT* matrix2)
{
  bool result = true;

  assert(matrix1 != NULL);
  assert(matrix2 != NULL);

  if (matrix1->numRows != matrix2->numColumns)
    return false;
  if (matrix1->numColumns != matrix2->numRows)
    return false;
  if (matrix1->numNonzeros != matrix2->numNonzeros)
    return false;

  int* currentColumnEntries = (int*) malloc(matrix1->numColumns * sizeof(int) );
  for (int column = 0; column < matrix2->numRows; ++column)
    currentColumnEntries[column] = matrix2->rowStarts[column];

  for (int row = 0; row < matrix1->numRows; ++row)
  {
    int begin = matrix1->rowStarts[row];
    int end = row + 1 < matrix1->numRows ? matrix1->rowStarts[row + 1] : matrix1->numNonzeros;
    for (int entry1 = begin; entry1 < end; ++entry1)
    {
      int column = matrix1->entryColumns[entry1];
      int entry2 = currentColumnEntries[column];
      if (matrix2->entryColumns[entry2] != row
        || matrix2->entryValues[entry2] != matrix1->entryValues[entry1])
      {
        result = false;
        goto cleanup;
      }
      currentColumnEntries[column]++;
    }
  }

cleanup:

  free(currentColumnEntries);

  return result;
}

bool TUchrmatCheckTranspose(TU_CHRMAT* matrix1, TU_CHRMAT* matrix2)
{
  bool result = true;

  assert(matrix1 != NULL);
  assert(matrix2 != NULL);

  if (matrix1->numRows != matrix2->numColumns)
    return false;
  if (matrix1->numColumns != matrix2->numRows)
    return false;
  if (matrix1->numNonzeros != matrix2->numNonzeros)
    return false;

  int* currentColumnEntries = (int*) malloc(matrix1->numColumns * sizeof(int) );
  for (int column = 0; column < matrix2->numRows; ++column)
    currentColumnEntries[column] = matrix2->rowStarts[column];

  for (int row = 0; row < matrix1->numRows; ++row)
  {
    int begin = matrix1->rowStarts[row];
    int end = row + 1 < matrix1->numRows ? matrix1->rowStarts[row + 1] : matrix1->numNonzeros;
    for (int entry1 = begin; entry1 < end; ++entry1)
    {
      int column = matrix1->entryColumns[entry1];
      int entry2 = currentColumnEntries[column];
      if (matrix2->entryColumns[entry2] != row
        || matrix2->entryValues[entry2] != matrix1->entryValues[entry1])
      {
        result = false;
        goto cleanup;
      }
      currentColumnEntries[column]++;
    }
  }

cleanup:

  free(currentColumnEntries);

  return result;
}

bool TUdblmatCheckSorted(TU_DBLMAT* sparse)
{
  assert(sparse != NULL);

  for (int row = 0; row < sparse->numRows; ++row)
  {
    int start = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row+1] : sparse->numNonzeros;
    for (int i = start + 1; i < end; ++i)
    {
      if (sparse->entryColumns[i-1] > sparse->entryColumns[i])
        return false;
    }
  }

  return true;
}

bool TUintmatCheckSorted(TU_INTMAT* sparse)
{
  return TUdblmatCheckSorted((TU_DBLMAT*) sparse);
}

bool TUchrmatCheckSorted(TU_CHRMAT* sparse)
{
  return TUdblmatCheckSorted((TU_DBLMAT*) sparse);
}

bool TUisBinaryDbl(TU* tu, TU_DBLMAT* sparse, double epsilon, TU_SUBMAT** submatrix)
{
  assert(sparse != NULL);

  for (int row = 0; row < sparse->numRows; ++row)
  {
    int begin = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      double value = sparse->entryValues[entry];
      int rounded = (int)(value + 0.5);
      if (rounded < 0 || rounded > +1 || fabs(value - rounded) > epsilon)
      {
        if (submatrix)
          TUsubmatCreate1x1(tu, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool TUisBinaryInt(TU* tu, TU_INTMAT* sparse, TU_SUBMAT** submatrix)
{
  assert(sparse != NULL);

  for (int row = 0; row < sparse->numRows; ++row)
  {
    int begin = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      int value = sparse->entryValues[entry];
      if (value < 0 || value > 1)
      {
        if (submatrix)
          TUsubmatCreate1x1(tu, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool TUisBinaryChr(TU* tu, TU_CHRMAT* sparse, TU_SUBMAT** submatrix)
{
  assert(sparse != NULL);

  for (int row = 0; row < sparse->numRows; ++row)
  {
    int begin = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      char value = sparse->entryValues[entry];
      if (value < 0 || value > 1)
      {
        if (submatrix)
          TUsubmatCreate1x1(tu, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool TUisTernaryDbl(TU* tu, TU_DBLMAT* sparse, double epsilon, TU_SUBMAT** submatrix)
{
  assert(sparse != NULL);

  for (int row = 0; row < sparse->numRows; ++row)
  {
    int begin = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      double value = sparse->entryValues[entry];
      int rounded = (int)(value + 0.5);
      if (rounded < -1 || rounded > +1 || fabs(value - rounded) > epsilon)
      {
        if (submatrix)
          TUsubmatCreate1x1(tu, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool TUisTernaryInt(TU* tu, TU_INTMAT* sparse, TU_SUBMAT** submatrix)
{
  assert(sparse != NULL);

  for (int row = 0; row < sparse->numRows; ++row)
  {
    int begin = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      int value = sparse->entryValues[entry];
      if (value < -1 || value > +1)
      {
        if (submatrix)
          TUsubmatCreate1x1(tu, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool TUisTernaryChr(TU* tu, TU_CHRMAT* sparse, TU_SUBMAT** submatrix)
{
  assert(sparse != NULL);

  for (int row = 0; row < sparse->numRows; ++row)
  {
    int begin = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      char value = sparse->entryValues[entry];
      if (value < -1 || value > +1)
      {
        if (submatrix)
          TUsubmatCreate1x1(tu, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

TU_ERROR TUsubmatCreate(TU* tu, TU_SUBMAT** psubmatrix, int numRows, int numColumns)
{
  assert(psubmatrix != NULL);

  TU_CALL( TUallocBlock(tu, psubmatrix) );
  TU_SUBMAT* submatrix = *psubmatrix;
  submatrix->numRows = numRows;
  submatrix->numColumns = numColumns;
  submatrix->rows = NULL;
  submatrix->columns = NULL;
  TU_CALL( TUallocBlockArray(tu, &submatrix->rows, numRows) );
  TU_CALL( TUallocBlockArray(tu, &submatrix->columns, numColumns) );

  return TU_OKAY;
}

void TUsubmatCreate1x1(TU* tu, TU_SUBMAT** submatrix, int row, int column)
{
  TUsubmatCreate(tu, submatrix, 1, 1);
  (*submatrix)->rows[0] = row;
  (*submatrix)->columns[0] = column;
}

TU_ERROR TUsubmatFree(TU* tu, TU_SUBMAT** psubmatrix)
{
  assert(psubmatrix);

  if ((*psubmatrix)->rows)
    TUfreeBlockArray(tu, &(*psubmatrix)->rows);
  if ((*psubmatrix)->columns)
    TUfreeBlockArray(tu, &(*psubmatrix)->columns);
  TUfreeBlockArray(tu, psubmatrix);
  *psubmatrix = NULL;

  return TU_OKAY;
}

static int TUsortSubmatrixCompare(const void* p1, const void* p2)
{
  return *(int*)p1 - *(int*)p2;
}

TU_ERROR TUsortSubmatrix(TU_SUBMAT* submatrix)
{
  assert(submatrix);

  qsort(submatrix->rows, submatrix->numRows, sizeof(int), TUsortSubmatrixCompare);
  qsort(submatrix->columns, submatrix->numColumns, sizeof(int), TUsortSubmatrixCompare);

  return TU_OKAY;
}

TU_ERROR TUdblmatFilterSubmat(TU* tu, TU_DBLMAT* matrix, TU_SUBMAT* submatrix, TU_DBLMAT** result)
{
  assert(matrix);
  assert(submatrix);
  assert(result);

  int* columnMap = NULL;
  TU_CALL( TUallocStackArray(tu, &columnMap, matrix->numColumns) );
  for (int c = 0; c < matrix->numColumns; ++c)
    columnMap[c] = -1;
  for (int j = 0; j < submatrix->numColumns; ++j)
  {
    assert(submatrix->columns[j] < matrix->numColumns);
    columnMap[submatrix->columns[j]] = j;
  }

  /* Count nonzeros. */
  int numNonzeros = 0;
  for (int i = 0; i < submatrix->numRows; ++i)
  {
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    int begin = matrix->rowStarts[r];
    int end = r + 1 < matrix->numRows ? matrix->rowStarts[r+1] : matrix->numNonzeros;
    for (int e = begin; e < end; ++e)
    {
      int c = matrix->entryColumns[e];
      if (columnMap[c] >= 0)
        ++numNonzeros;
    }
  }

  TU_CALL( TUdblmatCreate(tu, result, submatrix->numRows, submatrix->numColumns, numNonzeros) );

  /* Copy nonzeros. */
  (*result)->numNonzeros = 0;
  for (int i = 0; i < submatrix->numRows; ++i)
  {
    (*result)->rowStarts[i] = (*result)->numNonzeros;
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    int begin = matrix->rowStarts[r];
    int end = r + 1 < matrix->numRows ? matrix->rowStarts[r+1] : matrix->numNonzeros;
    for (int e = begin; e < end; ++e)
    {
      int c = matrix->entryColumns[e];
      if (columnMap[c] >= 0)
      {
        (*result)->entryColumns[(*result)->numNonzeros] = columnMap[c];
        (*result)->entryValues[(*result)->numNonzeros] = matrix->entryValues[e];
        (*result)->numNonzeros++;
      }
    }
  }
  (*result)->rowStarts[(*result)->numRows] = (*result)->numNonzeros;

  if (columnMap)
    TU_CALL( TUfreeStackArray(tu, &columnMap) );

  return TU_OKAY;
}


TU_ERROR TUintmatFilterSubmat(TU* tu, TU_INTMAT* matrix, TU_SUBMAT* submatrix, TU_INTMAT** result)
{
  assert(matrix);
  assert(submatrix);
  assert(result);

  int* columnMap = NULL;
  TU_CALL( TUallocStackArray(tu, &columnMap, matrix->numColumns) );
  for (int c = 0; c < matrix->numColumns; ++c)
    columnMap[c] = -1;
  for (int j = 0; j < submatrix->numColumns; ++j)
  {
    assert(submatrix->columns[j] < matrix->numColumns);
    columnMap[submatrix->columns[j]] = j;
  }

  /* Count nonzeros. */
  int numNonzeros = 0;
  for (int i = 0; i < submatrix->numRows; ++i)
  {
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    int begin = matrix->rowStarts[r];
    int end = r + 1 < matrix->numRows ? matrix->rowStarts[r+1] : matrix->numNonzeros;
    for (int e = begin; e < end; ++e)
    {
      int c = matrix->entryColumns[e];
      if (columnMap[c] >= 0)
        ++numNonzeros;
    }
  }

  TU_CALL( TUintmatCreate(tu, result, submatrix->numRows, submatrix->numColumns, numNonzeros) );

  /* Copy nonzeros. */
  (*result)->numNonzeros = 0;
  for (int i = 0; i < submatrix->numRows; ++i)
  {
    (*result)->rowStarts[i] = (*result)->numNonzeros;
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    int begin = matrix->rowStarts[r];
    int end = r + 1 < matrix->numRows ? matrix->rowStarts[r+1] : matrix->numNonzeros;
    for (int e = begin; e < end; ++e)
    {
      int c = matrix->entryColumns[e];
      if (columnMap[c] >= 0)
      {
        (*result)->entryColumns[(*result)->numNonzeros] = columnMap[c];
        (*result)->entryValues[(*result)->numNonzeros] = matrix->entryValues[e];
        (*result)->numNonzeros++;
      }
    }
  }
  (*result)->rowStarts[(*result)->numRows] = (*result)->numNonzeros;

  if (columnMap)
    TU_CALL( TUfreeStackArray(tu, &columnMap) );

  return TU_OKAY;
}

TU_ERROR TUchrmatFilterSubmat(TU* tu, TU_CHRMAT* matrix, TU_SUBMAT* submatrix, TU_CHRMAT** result)
{
  assert(matrix);
  assert(submatrix);
  assert(result);

  int* columnMap = NULL;
  TU_CALL( TUallocStackArray(tu, &columnMap, matrix->numColumns) );
  for (int c = 0; c < matrix->numColumns; ++c)
    columnMap[c] = -1;
  for (int j = 0; j < submatrix->numColumns; ++j)
  {
    assert(submatrix->columns[j] < matrix->numColumns);
    columnMap[submatrix->columns[j]] = j;
  }

  /* Count nonzeros. */
  int numNonzeros = 0;
  for (int i = 0; i < submatrix->numRows; ++i)
  {
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    int begin = matrix->rowStarts[r];
    int end = r + 1 < matrix->numRows ? matrix->rowStarts[r+1] : matrix->numNonzeros;
    for (int e = begin; e < end; ++e)
    {
      int c = matrix->entryColumns[e];
      if (columnMap[c] >= 0)
        ++numNonzeros;
    }
  }

  TU_CALL( TUchrmatCreate(tu, result, submatrix->numRows, submatrix->numColumns, numNonzeros) );

  /* Copy nonzeros. */
  (*result)->numNonzeros = 0;
  for (int i = 0; i < submatrix->numRows; ++i)
  {
    (*result)->rowStarts[i] = (*result)->numNonzeros;
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    int begin = matrix->rowStarts[r];
    int end = r + 1 < matrix->numRows ? matrix->rowStarts[r+1] : matrix->numNonzeros;
    for (int e = begin; e < end; ++e)
    {
      int c = matrix->entryColumns[e];
      if (columnMap[c] >= 0)
      {
        (*result)->entryColumns[(*result)->numNonzeros] = columnMap[c];
        (*result)->entryValues[(*result)->numNonzeros] = matrix->entryValues[e];
        (*result)->numNonzeros++;
      }
    }
  }
  (*result)->rowStarts[(*result)->numRows] = (*result)->numNonzeros;

  if (columnMap)
    TU_CALL( TUfreeStackArray(tu, &columnMap) );

  return TU_OKAY;
}

TU_ERROR TUsupportDbl(TU* tu, TU_DBLMAT* matrix, double epsilon, TU_CHRMAT** psupport)
{
  assert(tu);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  TU_CALL( TUchrmatCreate(tu, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  TU_CHRMAT* result = *psupport;

  int resultEntry = 0;
  for (int row = 0; row < matrix->numRows; ++row)
  {
    result->rowStarts[row] = resultEntry;
    int firstMatrixEntry = matrix->rowStarts[row];
    int beyondMatrixEntry = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int matrixEntry = firstMatrixEntry; matrixEntry < beyondMatrixEntry; ++matrixEntry)
    {
      if (fabs(matrix->entryValues[matrixEntry]) > epsilon)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = 1;
        resultEntry++;
      }
    }
  }
  result->rowStarts[matrix->numRows] = resultEntry;

  return TU_OKAY;
}

TU_ERROR TUsignedSupportDbl(TU* tu, TU_DBLMAT* matrix, double epsilon, TU_CHRMAT** psupport)
{
  assert(tu);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  TU_CALL( TUchrmatCreate(tu, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  TU_CHRMAT* result = *psupport;

  int resultEntry = 0;
  for (int row = 0; row < matrix->numRows; ++row)
  {
    result->rowStarts[row] = resultEntry;
    int firstMatrixEntry = matrix->rowStarts[row];
    int beyondMatrixEntry = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int matrixEntry = firstMatrixEntry; matrixEntry < beyondMatrixEntry; ++matrixEntry)
    {
      if (matrix->entryValues[matrixEntry] > epsilon)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = 1;
        resultEntry++;
      }
      else if (matrix->entryValues[matrixEntry] < -epsilon)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = -1;
        resultEntry++;
      }
    }
  }
  result->rowStarts[matrix->numRows] = resultEntry;

  return TU_OKAY;
}

TU_ERROR TUsupportInt(TU* tu, TU_INTMAT* matrix, TU_CHRMAT** psupport)
{
  assert(tu);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  TU_CALL( TUchrmatCreate(tu, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  TU_CHRMAT* result = *psupport;

  int resultEntry = 0;
  for (int row = 0; row < matrix->numRows; ++row)
  {
    result->rowStarts[row] = resultEntry;
    int firstMatrixEntry = matrix->rowStarts[row];
    int beyondMatrixEntry = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int matrixEntry = firstMatrixEntry; matrixEntry < beyondMatrixEntry; ++matrixEntry)
    {
      if (matrix->entryValues[matrixEntry] != 0)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = 1;
        resultEntry++;
      }
    }
  }
  result->rowStarts[matrix->numRows] = resultEntry;

  return TU_OKAY;
}

TU_ERROR TUsignedSupportInt(TU* tu, TU_INTMAT* matrix, TU_CHRMAT** psupport)
{
  assert(tu);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  TU_CALL( TUchrmatCreate(tu, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  TU_CHRMAT* result = *psupport;

  int resultEntry = 0;
  for (int row = 0; row < matrix->numRows; ++row)
  {
    result->rowStarts[row] = resultEntry;
    int firstMatrixEntry = matrix->rowStarts[row];
    int beyondMatrixEntry = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int matrixEntry = firstMatrixEntry; matrixEntry < beyondMatrixEntry; ++matrixEntry)
    {
      if (matrix->entryValues[matrixEntry] > 0)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = 1;
        resultEntry++;
      }
      else if (matrix->entryValues[matrixEntry] < 0)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = -1;
        resultEntry++;
      }
    }
  }
  result->rowStarts[matrix->numRows] = resultEntry;

  return TU_OKAY;
}

TU_ERROR TUsupportChr(TU* tu, TU_CHRMAT* matrix, TU_CHRMAT** psupport)
{
  assert(tu);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  TU_CALL( TUchrmatCreate(tu, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  TU_CHRMAT* result = *psupport;

  int resultEntry = 0;
  for (int row = 0; row < matrix->numRows; ++row)
  {
    result->rowStarts[row] = resultEntry;
    int firstMatrixEntry = matrix->rowStarts[row];
    int beyondMatrixEntry = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int matrixEntry = firstMatrixEntry; matrixEntry < beyondMatrixEntry; ++matrixEntry)
    {
      if (matrix->entryValues[matrixEntry] != 0)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = 1;
        resultEntry++;
      }
    }
  }
  result->rowStarts[matrix->numRows] = resultEntry;

  return TU_OKAY;
}

TU_ERROR TUsignedSupportChr(TU* tu, TU_CHRMAT* matrix, TU_CHRMAT** psupport)
{
  assert(tu);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  TU_CALL( TUchrmatCreate(tu, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  TU_CHRMAT* result = *psupport;

  int resultEntry = 0;
  for (int row = 0; row < matrix->numRows; ++row)
  {
    result->rowStarts[row] = resultEntry;
    int firstMatrixEntry = matrix->rowStarts[row];
    int beyondMatrixEntry = (row+1 < matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int matrixEntry = firstMatrixEntry; matrixEntry < beyondMatrixEntry; ++matrixEntry)
    {
      if (matrix->entryValues[matrixEntry] > 0)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = 1;
        resultEntry++;
      }
      else if (matrix->entryValues[matrixEntry] < 0)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = -1;
        resultEntry++;
      }
    }
  }
  result->rowStarts[matrix->numRows] = resultEntry;

  return TU_OKAY;
}

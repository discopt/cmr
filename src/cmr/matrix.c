#include <cmr/matrix.h>

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "sort.h"
#include "env_internal.h"

CMR_ERROR CMRdblmatCreate(CMR* cmr, CMR_DBLMAT** matrix, int numRows, int numColumns,
  int numNonzeros)
{
  assert(matrix);
  assert(*matrix == NULL);

  CMR_CALL( CMRallocBlock(cmr, matrix) );
  (*matrix)->numRows = numRows;
  (*matrix)->numColumns = numColumns;
  (*matrix)->numNonzeros = numNonzeros;
  (*matrix)->rowStarts = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->rowStarts, numRows + 1) );
  if (numNonzeros > 0)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryColumns, numNonzeros) );
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryValues, numNonzeros) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatFree(CMR* cmr, CMR_DBLMAT** matrix)
{
  assert(matrix);
  assert(*matrix);
  assert((*matrix)->rowStarts);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryColumns);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryValues);

  CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->rowStarts) );
  if ((*matrix)->numNonzeros > 0)
  {
    CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->entryColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->entryValues) );
  }
  CMR_CALL( CMRfreeBlock(cmr, matrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatChangeNumNonzeros(CMR* cmr, CMR_DBLMAT* matrix, int newNumNonzeros)
{
  assert(cmr);
  assert(matrix);

  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryColumns, newNumNonzeros) );
  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatCopy(CMR* cmr, CMR_DBLMAT* matrix, CMR_DBLMAT** result)
{
  assert(cmr);
  assert(matrix);
  assert(result);
  assert(*result == NULL);

  CMR_CALL( CMRdblmatCreate(cmr, result, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  for (int row = 0; row < matrix->numRows; ++row)
    (*result)->rowStarts[row] = matrix->rowStarts[row];
  (*result)->rowStarts[matrix->numRows] = matrix->numNonzeros;
  for (int entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    (*result)->entryColumns[entry] = matrix->entryColumns[entry];
    (*result)->entryValues[entry] = matrix->entryValues[entry];
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatTranspose(CMR* cmr, CMR_DBLMAT* matrix, CMR_DBLMAT** result)
{
  assert(cmr);
  assert(matrix);
  assert(result);
  assert(*result == NULL);
  assert(CMRdblmatCheckSorted(matrix));

  CMR_CALL( CMRdblmatCreate(cmr, result, matrix->numColumns, matrix->numRows, matrix->numNonzeros) );

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

  return CMR_OKAY;
}

CMR_ERROR CMRintmatCreate(CMR* cmr, CMR_INTMAT** matrix, int numRows, int numColumns, int numNonzeros)
{
  assert(matrix);
  assert(*matrix == NULL);

  CMR_CALL( CMRallocBlock(cmr, matrix) );
  (*matrix)->numRows = numRows;
  (*matrix)->numColumns = numColumns;
  (*matrix)->numNonzeros = numNonzeros;
  (*matrix)->rowStarts = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->rowStarts, numRows + 1) );
  if (numNonzeros > 0)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryColumns, numNonzeros) );
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryValues, numNonzeros) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRintmatFree(CMR* cmr, CMR_INTMAT** matrix)
{
  assert(matrix);
  assert(*matrix);
  assert((*matrix)->rowStarts);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryColumns);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryValues);

  CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->rowStarts) );
  if ((*matrix)->numNonzeros > 0)
  {
    CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->entryColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->entryValues) );
  }
  CMR_CALL( CMRfreeBlock(cmr, matrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRintmatChangeNumNonzeros(CMR* cmr, CMR_INTMAT* matrix, int newNumNonzeros)
{
  assert(cmr);
  assert(matrix);

  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryColumns, newNumNonzeros) );
  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return CMR_OKAY;
}

CMR_ERROR CMRintmatCopy(CMR* cmr, CMR_INTMAT* matrix, CMR_INTMAT** result)
{
  assert(cmr);
  assert(matrix);
  assert(result);
  assert(*result == NULL);

  CMR_CALL( CMRintmatCreate(cmr, result, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  for (int row = 0; row < matrix->numRows; ++row)
    (*result)->rowStarts[row] = matrix->rowStarts[row];
  (*result)->rowStarts[matrix->numRows] = matrix->numNonzeros;
  for (int entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    (*result)->entryColumns[entry] = matrix->entryColumns[entry];
    (*result)->entryValues[entry] = matrix->entryValues[entry];
  }

  return CMR_OKAY;
}

CMR_ERROR CMRintmatTranspose(CMR* cmr, CMR_INTMAT* matrix, CMR_INTMAT** result)
{
  assert(cmr);
  assert(matrix);
  assert(result);
  assert(*result == NULL);
  assert(CMRintmatCheckSorted(matrix));

  CMR_CALL( CMRintmatCreate(cmr, result, matrix->numColumns, matrix->numRows, matrix->numNonzeros) );

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

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatCreate(CMR* cmr, CMR_CHRMAT** matrix, int numRows, int numColumns, int numNonzeros)
{
  assert(matrix);
  assert(*matrix == NULL);

  CMR_CALL( CMRallocBlock(cmr, matrix) );
  (*matrix)->numRows = numRows;
  (*matrix)->numColumns = numColumns;
  (*matrix)->numNonzeros = numNonzeros;
  (*matrix)->rowStarts = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->rowStarts, numRows + 1) );
  if (numNonzeros > 0)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryColumns, numNonzeros) );
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryValues, numNonzeros) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatFree(CMR* cmr, CMR_CHRMAT** matrix)
{
  assert(matrix);
  assert(*matrix);
  assert((*matrix)->rowStarts);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryColumns);
  assert((*matrix)->numNonzeros == 0 || (*matrix)->entryValues);

  CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->rowStarts) );
  if ((*matrix)->numNonzeros > 0)
  {
    CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->entryColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &(*matrix)->entryValues) );
  }
  CMR_CALL( CMRfreeBlock(cmr, matrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatChangeNumNonzeros(CMR* cmr, CMR_CHRMAT* matrix, int newNumNonzeros)
{
  assert(cmr);
  assert(matrix);

  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryColumns, newNumNonzeros) );
  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatCopy(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT** result)
{
  assert(cmr);
  assert(matrix);
  assert(result);
  assert(*result == NULL);

  CMR_CALL( CMRchrmatCreate(cmr, result, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  for (int row = 0; row < matrix->numRows; ++row)
    (*result)->rowStarts[row] = matrix->rowStarts[row];
  (*result)->rowStarts[matrix->numRows] = matrix->numNonzeros;
  for (int entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    (*result)->entryColumns[entry] = matrix->entryColumns[entry];
    (*result)->entryValues[entry] = matrix->entryValues[entry];
  }

  return CMR_OKAY;
}


CMR_ERROR CMRchrmatTranspose(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT** result)
{
  assert(cmr);
  assert(matrix);
  assert(result);
  assert(*result == NULL);
  assert(CMRchrmatCheckSorted(matrix));

  CMR_CALL( CMRchrmatCreate(cmr, result, matrix->numColumns, matrix->numRows, matrix->numNonzeros) );

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

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatPrintSparse(FILE* stream, CMR_DBLMAT* matrix)
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

  return CMR_OKAY;
}

CMR_ERROR CMRintmatPrintSparse(FILE* stream, CMR_INTMAT* matrix)
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

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatPrintSparse(FILE* stream, CMR_CHRMAT* matrix)
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

  return CMR_OKAY;
}



CMR_ERROR CMRdblmatPrintDense(FILE* stream, CMR_DBLMAT* matrix, char zeroChar, bool header)
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

  return CMR_OKAY;
}

CMR_ERROR CMRintmatPrintDense(FILE* stream, CMR_INTMAT* matrix, char zeroChar, bool header)
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

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatPrintDense(FILE* stream, CMR_CHRMAT* matrix, char zeroChar, bool header)
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

  return CMR_OKAY;
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

CMR_ERROR CMRdblmatCreateFromSparseStream(CMR* cmr, CMR_DBLMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%d %d %d", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return CMR_ERROR_INPUT;

  /* Read all nonzeros. */

  DblNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
  int entry = 0;
  for (int i = 0; i < numNonzeros; ++i)
  {
    int row;
    int column;
    double value;
    numRead = fscanf(stream, "%d %d %lf", &row, &column, &value);
    if (numRead < 3 || row < 0 || column < 0 || row >= numRows || column >= numColumns)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      return CMR_ERROR_INPUT;
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
  CMR_CALL( CMRsort(cmr, numNonzeros, nonzeros, sizeof(DblNonzero), compareDblNonzeros) );

  CMR_CALL( CMRdblmatCreate(cmr, pmatrix, numRows, numColumns, numNonzeros) );
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
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      CMR_CALL( CMRdblmatFree(cmr, pmatrix) );
      return CMR_ERROR_INPUT;
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

  CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );

  return CMR_OKAY;
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

CMR_ERROR CMRintmatCreateFromSparseStream(CMR* cmr, CMR_INTMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%d %d %d", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return CMR_ERROR_INPUT;

  /* Read all nonzeros. */

  IntNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
  int entry = 0;
  for (int i = 0; i < numNonzeros; ++i)
  {
    int row;
    int column;
    int value;
    numRead = fscanf(stream, "%d %d %d", &row, &column, &value);
    if (numRead < 3 || row < 0 || column < 0 || row >= numRows || column >= numColumns)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      return CMR_ERROR_INPUT;
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
  CMR_CALL( CMRsort(cmr, numNonzeros, nonzeros, sizeof(IntNonzero), compareIntNonzeros) );

  CMR_CALL( CMRintmatCreate(cmr, pmatrix, numRows, numColumns, numNonzeros) );
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
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      CMR_CALL( CMRintmatFree(cmr, pmatrix) );
      return CMR_ERROR_INPUT;
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

  CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );

  return CMR_OKAY;
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

CMR_ERROR CMRchrmatCreateFromSparseStream(CMR* cmr, CMR_CHRMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%d %d %d", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return CMR_ERROR_INPUT;

  /* Read all nonzeros. */

  ChrNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
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
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      return CMR_ERROR_INPUT;
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
  CMR_CALL( CMRsort(cmr, numNonzeros, nonzeros, sizeof(ChrNonzero), compareChrNonzeros) );

  CMR_CALL( CMRchrmatCreate(cmr, pmatrix, numRows, numColumns, numNonzeros) );
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
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      CMR_CALL( CMRchrmatFree(cmr, pmatrix) );
      return CMR_ERROR_INPUT;
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

  CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatCreateFromDenseStream(CMR* cmr, CMR_DBLMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns;
  int numRead = fscanf(stream, "%d %d", &numRows, &numColumns);
  if (numRead < 2)
    return CMR_ERROR_INPUT;

  CMR_CALL( CMRdblmatCreate(cmr, pmatrix, numRows, numColumns, 0) );

  /* Initial memory. */
  int memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  int* entryColumns = NULL;
  double* entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &entryColumns, memEntries) );
  CMR_CALL( CMRallocBlockArray(cmr, &entryValues, memEntries) );

  int entry = 0;
  for (int row = 0; row < numRows; ++row)
  {
    (*pmatrix)->rowStarts[row] = entry;
    for (int column = 0; column < numColumns; ++column)
    {
      double x;
      numRead = fscanf(stream, "%lf", &x);
      if (numRead < 1)
        return CMR_ERROR_INPUT;

      if (x == 0.0)
        continue;

      if (entry == memEntries)
      {
        memEntries = 2*memEntries;
        CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, memEntries) );
        CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, memEntries) );
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
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry) );
  }

  (*pmatrix)->entryColumns = entryColumns;
  (*pmatrix)->entryValues = entryValues;
  (*pmatrix)->numNonzeros = entry;

  return CMR_OKAY;
}

CMR_ERROR CMRintmatCreateFromDenseStream(CMR* cmr, CMR_INTMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows, numColumns;
  int numRead = fscanf(stream, "%d %d", &numRows, &numColumns);
  if (numRead < 2)
    return CMR_ERROR_INPUT;

  CMR_CALL( CMRintmatCreate(cmr, pmatrix, numRows, numColumns, 0) );

  /* Initial memory. */
  int memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  int* entryColumns = NULL;
  int* entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &entryColumns, memEntries) );
  CMR_CALL( CMRallocBlockArray(cmr, &entryValues, memEntries) );

  int entry = 0;
  for (int row = 0; row < numRows; ++row)
  {
    (*pmatrix)->rowStarts[row] = entry;
    for (int column = 0; column < numColumns; ++column)
    {
      int x;
      numRead = fscanf(stream, "%d", &x);
      if (numRead < 1)
        return CMR_ERROR_INPUT;

      if (x == 0.0)
        continue;

      if (entry == memEntries)
      {
        memEntries = 2*memEntries;
        CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, memEntries) );
        CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, memEntries) );
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
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry) );
  }

  (*pmatrix)->entryColumns = entryColumns;
  (*pmatrix)->entryValues = entryValues;
  (*pmatrix)->numNonzeros = entry;

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatCreateFromDenseStream(CMR* cmr, CMR_CHRMAT** pmatrix, FILE* stream)
{
  assert(pmatrix);
  assert(!*pmatrix);
  assert(stream);

  int numRows = 0;
  int numColumns = 0;
  int numRead = fscanf(stream, "%d %d", &numRows, &numColumns);
  if (numRead < 2)
    return CMR_ERROR_INPUT;

  CMR_CALL( CMRchrmatCreate(cmr, pmatrix, numRows, numColumns, 0) );

  /* Initial memory. */
  int memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  int* entryColumns = NULL;
  char* entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &entryColumns, memEntries) );
  CMR_CALL( CMRallocBlockArray(cmr, &entryValues, memEntries) );

  int entry = 0;
  for (int row = 0; row < numRows; ++row)
  {
    (*pmatrix)->rowStarts[row] = entry;
    for (int column = 0; column < numColumns; ++column)
    {
      int x;
      numRead = fscanf(stream, "%d", &x);
      if (numRead < 1)
        return CMR_ERROR_INPUT;

      if (x == 0.0)
        continue;

      if (entry == memEntries)
      {
        memEntries = 2*memEntries;
        CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, memEntries) );
        CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, memEntries) );
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
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry) );
  }

  (*pmatrix)->entryColumns = entryColumns;
  (*pmatrix)->entryValues = entryValues;
  (*pmatrix)->numNonzeros = entry;

  return CMR_OKAY;
}

bool CMRdblmatCheckEqual(CMR_DBLMAT* matrix1, CMR_DBLMAT* matrix2)
{
  assert(CMRdblmatCheckSorted(matrix1));
  assert(CMRdblmatCheckSorted(matrix2));

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

bool CMRintmatCheckEqual(CMR_INTMAT* matrix1, CMR_INTMAT* matrix2)
{
  assert(CMRintmatCheckSorted(matrix1));
  assert(CMRintmatCheckSorted(matrix2));

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

bool CMRchrmatCheckEqual(CMR_CHRMAT* matrix1, CMR_CHRMAT* matrix2)
{
  assert(CMRchrmatCheckSorted(matrix1));
  assert(CMRchrmatCheckSorted(matrix2));

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

bool CMRdblmatCheckTranspose(CMR_DBLMAT* matrix1, CMR_DBLMAT* matrix2)
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

bool CMRintmatCheckTranspose(CMR_INTMAT* matrix1, CMR_INTMAT* matrix2)
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

bool CMRchrmatCheckTranspose(CMR_CHRMAT* matrix1, CMR_CHRMAT* matrix2)
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

bool CMRdblmatCheckSorted(CMR_DBLMAT* sparse)
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

bool CMRintmatCheckSorted(CMR_INTMAT* sparse)
{
  return CMRdblmatCheckSorted((CMR_DBLMAT*) sparse);
}

bool CMRchrmatCheckSorted(CMR_CHRMAT* sparse)
{
  return CMRdblmatCheckSorted((CMR_DBLMAT*) sparse);
}

bool CMRisBinaryDbl(CMR* cmr, CMR_DBLMAT* sparse, double epsilon, CMR_SUBMAT** submatrix)
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
          CMRsubmatCreate1x1(cmr, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool CMRisBinaryInt(CMR* cmr, CMR_INTMAT* sparse, CMR_SUBMAT** submatrix)
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
          CMRsubmatCreate1x1(cmr, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool CMRisBinaryChr(CMR* cmr, CMR_CHRMAT* sparse, CMR_SUBMAT** submatrix)
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
          CMRsubmatCreate1x1(cmr, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool CMRisTernaryDbl(CMR* cmr, CMR_DBLMAT* sparse, double epsilon, CMR_SUBMAT** submatrix)
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
          CMRsubmatCreate1x1(cmr, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool CMRisTernaryInt(CMR* cmr, CMR_INTMAT* sparse, CMR_SUBMAT** submatrix)
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
          CMRsubmatCreate1x1(cmr, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool CMRisTernaryChr(CMR* cmr, CMR_CHRMAT* sparse, CMR_SUBMAT** submatrix)
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
          CMRsubmatCreate1x1(cmr, submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

CMR_ERROR CMRsubmatCreate(CMR* cmr, CMR_SUBMAT** psubmatrix, int numRows, int numColumns)
{
  assert(psubmatrix != NULL);

  CMR_CALL( CMRallocBlock(cmr, psubmatrix) );
  CMR_SUBMAT* submatrix = *psubmatrix;
  submatrix->numRows = numRows;
  submatrix->numColumns = numColumns;
  submatrix->rows = NULL;
  submatrix->columns = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &submatrix->rows, numRows) );
  CMR_CALL( CMRallocBlockArray(cmr, &submatrix->columns, numColumns) );

  return CMR_OKAY;
}

void CMRsubmatCreate1x1(CMR* cmr, CMR_SUBMAT** submatrix, int row, int column)
{
  CMRsubmatCreate(cmr, submatrix, 1, 1);
  (*submatrix)->rows[0] = row;
  (*submatrix)->columns[0] = column;
}

CMR_ERROR CMRsubmatFree(CMR* cmr, CMR_SUBMAT** psubmatrix)
{
  assert(psubmatrix);

  if ((*psubmatrix)->rows)
    CMRfreeBlockArray(cmr, &(*psubmatrix)->rows);
  if ((*psubmatrix)->columns)
    CMRfreeBlockArray(cmr, &(*psubmatrix)->columns);
  CMRfreeBlockArray(cmr, psubmatrix);
  *psubmatrix = NULL;

  return CMR_OKAY;
}

static int CMRsortSubmatrixCompare(const void* p1, const void* p2)
{
  return *(int*)p1 - *(int*)p2;
}

CMR_ERROR CMRsortSubmatrix(CMR_SUBMAT* submatrix)
{
  assert(submatrix);

  qsort(submatrix->rows, submatrix->numRows, sizeof(int), CMRsortSubmatrixCompare);
  qsort(submatrix->columns, submatrix->numColumns, sizeof(int), CMRsortSubmatrixCompare);

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatFilterSubmat(CMR* cmr, CMR_DBLMAT* matrix, CMR_SUBMAT* submatrix, CMR_DBLMAT** result)
{
  assert(matrix);
  assert(submatrix);
  assert(result);

  int* columnMap = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnMap, matrix->numColumns) );
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

  CMR_CALL( CMRdblmatCreate(cmr, result, submatrix->numRows, submatrix->numColumns, numNonzeros) );

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
    CMR_CALL( CMRfreeStackArray(cmr, &columnMap) );

  return CMR_OKAY;
}


CMR_ERROR CMRintmatFilterSubmat(CMR* cmr, CMR_INTMAT* matrix, CMR_SUBMAT* submatrix, CMR_INTMAT** result)
{
  assert(matrix);
  assert(submatrix);
  assert(result);

  int* columnMap = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnMap, matrix->numColumns) );
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

  CMR_CALL( CMRintmatCreate(cmr, result, submatrix->numRows, submatrix->numColumns, numNonzeros) );

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
    CMR_CALL( CMRfreeStackArray(cmr, &columnMap) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatFilterSubmat(CMR* cmr, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix, CMR_CHRMAT** result)
{
  assert(matrix);
  assert(submatrix);
  assert(result);

  int* columnMap = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnMap, matrix->numColumns) );
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

  CMR_CALL( CMRchrmatCreate(cmr, result, submatrix->numRows, submatrix->numColumns, numNonzeros) );

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
    CMR_CALL( CMRfreeStackArray(cmr, &columnMap) );

  return CMR_OKAY;
}

CMR_ERROR CMRsupportDbl(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, CMR_CHRMAT** psupport)
{
  assert(cmr);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  CMR_CALL( CMRchrmatCreate(cmr, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *psupport;

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

  return CMR_OKAY;
}

CMR_ERROR CMRsignedSupportDbl(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, CMR_CHRMAT** psupport)
{
  assert(cmr);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  CMR_CALL( CMRchrmatCreate(cmr, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *psupport;

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

  return CMR_OKAY;
}

CMR_ERROR CMRsupportInt(CMR* cmr, CMR_INTMAT* matrix, CMR_CHRMAT** psupport)
{
  assert(cmr);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  CMR_CALL( CMRchrmatCreate(cmr, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *psupport;

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

  return CMR_OKAY;
}

CMR_ERROR CMRsignedSupportInt(CMR* cmr, CMR_INTMAT* matrix, CMR_CHRMAT** psupport)
{
  assert(cmr);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  CMR_CALL( CMRchrmatCreate(cmr, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *psupport;

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

  return CMR_OKAY;
}

CMR_ERROR CMRsupportChr(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT** psupport)
{
  assert(cmr);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  CMR_CALL( CMRchrmatCreate(cmr, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *psupport;

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

  return CMR_OKAY;
}

CMR_ERROR CMRsignedSupportChr(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT** psupport)
{
  assert(cmr);
  assert(matrix);
  assert(psupport);
  assert(!*psupport);

  CMR_CALL( CMRchrmatCreate(cmr, psupport, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *psupport;

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

  return CMR_OKAY;
}

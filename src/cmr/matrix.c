// #define CMR_DEBUG

#include <cmr/matrix.h>

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>

#include "sort.h"
#include "env_internal.h"

CMR_ERROR CMRsubmatCreate(CMR* cmr, size_t numRows, size_t numColumns, CMR_SUBMAT** psubmatrix)
{
  assert(psubmatrix);

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

CMR_ERROR CMRsubmatCreate1x1(CMR* cmr, size_t row, size_t column, CMR_SUBMAT** psubmatrix)
{
  assert(psubmatrix);

  CMR_CALL( CMRsubmatCreate(cmr, 1, 1, psubmatrix) );
  CMR_SUBMAT* submatrix = *psubmatrix;
  submatrix->rows[0] = row;
  submatrix->columns[0] = column;

  return CMR_OKAY;
}

CMR_ERROR CMRsubmatFree(CMR* cmr, CMR_SUBMAT** psubmatrix)
{
  assert(psubmatrix);
  CMR_SUBMAT* submatrix = *psubmatrix;
  if (!submatrix)
    return CMR_OKAY;

  if (submatrix->rows)
    CMRfreeBlockArray(cmr, &submatrix->rows);
  if (submatrix->columns)
    CMRfreeBlockArray(cmr, &submatrix->columns);
  CMRfreeBlockArray(cmr, psubmatrix);

  return CMR_OKAY;
}

static int CMRsortSubmatrixCompare(const void* p1, const void* p2)
{
  return *(size_t*)p1 - *(size_t*)p2;
}

CMR_ERROR CMRsortSubmatrix(CMR* cmr, CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(submatrix);

  CMR_CALL( CMRsort(cmr, submatrix->numRows, submatrix->rows, sizeof(size_t), CMRsortSubmatrixCompare) );
  CMR_CALL( CMRsort(cmr, submatrix->numColumns, submatrix->columns, sizeof(size_t), CMRsortSubmatrixCompare) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatCreate(CMR* cmr, CMR_DBLMAT** matrix, int numRows, int numColumns,
  int numNonzeros)
{
  assert(matrix);
  assert(*matrix == NULL);

  CMR_CALL( CMRallocBlock(cmr, matrix) );
  (*matrix)->numRows = numRows;
  (*matrix)->numColumns = numColumns;
  (*matrix)->numNonzeros = numNonzeros;
  (*matrix)->rowSlice = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->rowSlice, numRows + 1) );
  if (numNonzeros > 0)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryColumns, numNonzeros) );
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryValues, numNonzeros) );
  }

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
  (*matrix)->rowSlice = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->rowSlice, numRows + 1) );
  if (numNonzeros > 0)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryColumns, numNonzeros) );
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryValues, numNonzeros) );
  }

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
  (*matrix)->rowSlice = NULL;
  (*matrix)->entryColumns = NULL;
  (*matrix)->entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->rowSlice, numRows + 1) );
  if (numNonzeros > 0)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryColumns, numNonzeros) );
    CMR_CALL( CMRallocBlockArray(cmr, &(*matrix)->entryValues, numNonzeros) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatFree(CMR* cmr, CMR_DBLMAT** pmatrix)
{
  assert(pmatrix);
  CMR_DBLMAT* matrix = *pmatrix;
  if (!matrix)
    return CMR_OKAY;

  assert(matrix->rowSlice);
  assert(matrix->numNonzeros == 0 || matrix->entryColumns);
  assert(matrix->numNonzeros == 0 || matrix->entryValues);

  CMR_CALL( CMRfreeBlockArray(cmr, &matrix->rowSlice) );
  if (matrix->numNonzeros > 0)
  {
    CMR_CALL( CMRfreeBlockArray(cmr, &matrix->entryColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &matrix->entryValues) );
  }
  CMR_CALL( CMRfreeBlock(cmr, pmatrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRintmatFree(CMR* cmr, CMR_INTMAT** pmatrix)
{
  assert(pmatrix);
  CMR_INTMAT* matrix = *pmatrix;
  if (!matrix)
    return CMR_OKAY;

  assert(matrix->rowSlice);
  assert(matrix->numNonzeros == 0 || matrix->entryColumns);
  assert(matrix->numNonzeros == 0 || matrix->entryValues);

  CMR_CALL( CMRfreeBlockArray(cmr, &matrix->rowSlice) );
  if (matrix->numNonzeros > 0)
  {
    CMR_CALL( CMRfreeBlockArray(cmr, &matrix->entryColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &matrix->entryValues) );
  }
  CMR_CALL( CMRfreeBlock(cmr, pmatrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatFree(CMR* cmr, CMR_CHRMAT** pmatrix)
{
  assert(pmatrix);
  CMR_CHRMAT* matrix = *pmatrix;
  if (!matrix)
    return CMR_OKAY;

  assert(matrix->rowSlice);
  assert(matrix->numNonzeros == 0 || matrix->entryColumns);
  assert(matrix->numNonzeros == 0 || matrix->entryValues);

  CMR_CALL( CMRfreeBlockArray(cmr, &matrix->rowSlice) );
  if (matrix->numNonzeros > 0)
  {
    CMR_CALL( CMRfreeBlockArray(cmr, &matrix->entryColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &matrix->entryValues) );
  }
  CMR_CALL( CMRfreeBlock(cmr, pmatrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatChangeNumNonzeros(CMR* cmr, CMR_DBLMAT* matrix, size_t newNumNonzeros)
{
  assert(cmr);
  assert(matrix);

  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryColumns, newNumNonzeros) );
  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return CMR_OKAY;
}

CMR_ERROR CMRintmatChangeNumNonzeros(CMR* cmr, CMR_INTMAT* matrix, size_t newNumNonzeros)
{
  assert(cmr);
  assert(matrix);

  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryColumns, newNumNonzeros) );
  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatChangeNumNonzeros(CMR* cmr, CMR_CHRMAT* matrix, size_t newNumNonzeros)
{
  assert(cmr);
  assert(matrix);

  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryColumns, newNumNonzeros) );
  CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryValues, newNumNonzeros) );
  matrix->numNonzeros = newNumNonzeros;

  return CMR_OKAY;
}

static
int compareEntries(const void** pa, const void** pb)
{
  size_t** a = (size_t**)(pa);
  size_t** b = (size_t**)(pb);
  return (**a < **b) ? -1 : (**a > **b ? 1 : 0);
}

CMR_ERROR CMRdblmatSortNonzeros(CMR* cmr, CMR_DBLMAT* matrix)
{
  assert(cmr);
  CMRconsistencyAssert( CMRdblmatConsistency(matrix) );

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    CMR_CALL( CMRsort2(cmr, beyond - first, &matrix->entryColumns[first], sizeof(size_t), &matrix->entryValues[first],
      sizeof(double), compareEntries) );
  }

  return CMR_OKAY;
}


CMR_ERROR CMRintmatSortNonzeros(CMR* cmr, CMR_INTMAT* matrix)
{
  assert(cmr);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    CMR_CALL( CMRsort2(cmr, beyond - first, &matrix->entryColumns[first], sizeof(size_t), &matrix->entryValues[first],
      sizeof(int), compareEntries) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatSortNonzeros(CMR* cmr, CMR_CHRMAT* matrix)
{
  assert(cmr);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    CMR_CALL( CMRsort2(cmr, beyond - first, &matrix->entryColumns[first], sizeof(size_t), &matrix->entryValues[first],
      sizeof(char), compareEntries) );
  }

  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatCopy(CMR* cmr, CMR_DBLMAT* matrix, CMR_DBLMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(*presult == NULL);

  CMR_CALL( CMRdblmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_DBLMAT* result = *presult;
  for (size_t row = 0; row <= matrix->numRows; ++row)
    result->rowSlice[row] = matrix->rowSlice[row];
  for (size_t entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    result->entryColumns[entry] = matrix->entryColumns[entry];
    result->entryValues[entry] = matrix->entryValues[entry];
  }

  return CMR_OKAY;
}

CMR_ERROR CMRintmatCopy(CMR* cmr, CMR_INTMAT* matrix, CMR_INTMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(*presult == NULL);

  CMR_CALL( CMRintmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_INTMAT* result = *presult;
  for (size_t row = 0; row <= matrix->numRows; ++row)
    result->rowSlice[row] = matrix->rowSlice[row];
  for (size_t entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    result->entryColumns[entry] = matrix->entryColumns[entry];
    result->entryValues[entry] = matrix->entryValues[entry];
  }

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatCopy(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(*presult == NULL);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;
  for (size_t row = 0; row <= matrix->numRows; ++row)
    result->rowSlice[row] = matrix->rowSlice[row];
  for (size_t entry = 0; entry < matrix->numNonzeros; ++entry)
  {
    result->entryColumns[entry] = matrix->entryColumns[entry];
    result->entryValues[entry] = matrix->entryValues[entry];
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatTranspose(CMR* cmr, CMR_DBLMAT* matrix, CMR_DBLMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(*presult == NULL);
  CMRconsistencyAssert( CMRdblmatConsistency(matrix) );

  CMR_CALL( CMRdblmatCreate(cmr, presult, matrix->numColumns, matrix->numRows, matrix->numNonzeros) );
  CMR_DBLMAT* result = *presult;

  /* Count number of nonzeros in each column, storing in the next entry. */
  for (size_t c = 0; c <= matrix->numColumns; ++c)
    result->rowSlice[c] = 0;
  for (size_t e = 0; e < matrix->numNonzeros; ++e)
    result->rowSlice[matrix->entryColumns[e] + 1]++;

  /* Compute start indices for columns. */
  for (size_t c = 1; c < matrix->numColumns; ++c)
    result->rowSlice[c] += result->rowSlice[c-1];

  /* Create nonzeros. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t column = matrix->entryColumns[entry];
      size_t transEntry = result->rowSlice[column];
      result->entryColumns[transEntry] = row;
      result->entryValues[transEntry] = matrix->entryValues[entry];
      result->rowSlice[column]++;
    }
  }

  /* We shifted rowSlice of result, so we shift it back. */
  for (size_t c = matrix->numColumns; c > 0; --c)
    result->rowSlice[c] = result->rowSlice[c-1];
  result->rowSlice[0] = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRintmatTranspose(CMR* cmr, CMR_INTMAT* matrix, CMR_INTMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(*presult == NULL);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );

  CMR_CALL( CMRintmatCreate(cmr, presult, matrix->numColumns, matrix->numRows, matrix->numNonzeros) );
  CMR_INTMAT* result = *presult;

  /* Count number of nonzeros in each column, storing in the next entry. */
  for (size_t c = 0; c <= matrix->numColumns; ++c)
    result->rowSlice[c] = 0;
  for (size_t e = 0; e < matrix->numNonzeros; ++e)
    result->rowSlice[matrix->entryColumns[e] + 1]++;

  /* Compute start indices for columns. */
  for (size_t c = 1; c < matrix->numColumns; ++c)
    result->rowSlice[c] += result->rowSlice[c-1];

  /* Create nonzeros. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t column = matrix->entryColumns[entry];
      size_t transEntry = result->rowSlice[column];
      result->entryColumns[transEntry] = row;
      result->entryValues[transEntry] = matrix->entryValues[entry];
      result->rowSlice[column]++;
    }
  }

  /* We shifted rowSlice of result, so we shift it back. */
  for (size_t c = matrix->numColumns; c > 0; --c)
    result->rowSlice[c] = result->rowSlice[c-1];
  result->rowSlice[0] = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatTranspose(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(!*presult);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numColumns, matrix->numRows, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  /* Count number of nonzeros in each column, storing in the next entry. */
  for (size_t c = 0; c <= matrix->numColumns; ++c)
    result->rowSlice[c] = 0;
  for (size_t e = 0; e < matrix->numNonzeros; ++e)
    result->rowSlice[matrix->entryColumns[e] + 1]++;

  /* Compute start indices for columns. */
  for (size_t c = 1; c < matrix->numColumns; ++c)
    result->rowSlice[c] += result->rowSlice[c-1];

  /* Create nonzeros. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t column = matrix->entryColumns[entry];
      size_t transEntry = result->rowSlice[column];
      result->entryColumns[transEntry] = row;
      result->entryValues[transEntry] = matrix->entryValues[entry];
      result->rowSlice[column]++;
    }
  }

  /* We shifted rowSlice of result, so we shift it back. */
  for (size_t c = matrix->numColumns; c > 0; --c)
    result->rowSlice[c] = result->rowSlice[c-1];
  result->rowSlice[0] = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatPrintSparse(CMR* cmr, CMR_DBLMAT* matrix, FILE* stream)
{
  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%lu %lu %lu\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t entry = first; entry < beyond; ++entry)
      fprintf(stream, "%lu %lu  %f\n", row+1, matrix->entryColumns[entry]+1, matrix->entryValues[entry]);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRintmatPrintSparse(CMR* cmr, CMR_INTMAT* matrix, FILE* stream)
{
  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%lu %lu %lu\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t entry = first; entry < beyond; ++entry)
      fprintf(stream, "%lu %lu  %d\n", row+1, matrix->entryColumns[entry]+1, matrix->entryValues[entry]);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatPrintSparse(CMR* cmr, CMR_CHRMAT* matrix, FILE* stream)
{
  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%lu %lu %lu\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t entry = first; entry < beyond; ++entry)
      fprintf(stream, "%lu %lu  %d\n", row+1, matrix->entryColumns[entry]+1, matrix->entryValues[entry]);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatPrintDense(CMR* cmr, CMR_DBLMAT* matrix, FILE* stream, char zeroChar, bool header)
{
  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%lu %lu\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%lu ", (column+1) % 10);
    fputs("\n  ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%lu| ", (row+1) % 10);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    size_t column = 0;
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t entryColumn = matrix->entryColumns[entry];
      for (; column < entryColumn; ++column)
        fprintf(stream, "%c ", zeroChar);
      fprintf(stream, "%f ", matrix->entryValues[entry]);
      ++column;
    }
    for (; column < matrix->numColumns; ++column)
      fprintf(stream, "%c ", zeroChar);
    fputc('\n', stream);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRintmatPrintDense(CMR* cmr, CMR_INTMAT* matrix, FILE* stream, char zeroChar, bool header)
{
  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%lu %lu\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%lu ", (column+1) % 10);
    fputs("\n  ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%lu| ", (row+1) % 10);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    size_t column = 0;
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t entryColumn = matrix->entryColumns[entry];
      for (; column < entryColumn; ++column)
        fprintf(stream, "%c ", zeroChar);
      fprintf(stream, "%d ", matrix->entryValues[entry]);
      ++column;
    }
    for (; column < matrix->numColumns; ++column)
      fprintf(stream, "%c ", zeroChar);
    fputc('\n', stream);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatPrintDense(CMR* cmr, CMR_CHRMAT* matrix, FILE* stream, char zeroChar, bool header)
{
  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%lu %lu\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%lu ", (column+1) % 10);
    fputs("\n  ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%lu| ", (row+1) % 10);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    size_t column = 0;
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t entryColumn = matrix->entryColumns[entry];
      for (; column < entryColumn; ++column)
        fprintf(stream, "%c ", zeroChar);
      fprintf(stream, "%d ", matrix->entryValues[entry]);
      ++column;
    }
    for (; column < matrix->numColumns; ++column)
      fprintf(stream, "%c ", zeroChar);
    fputc('\n', stream);
  }

  return CMR_OKAY;
}



typedef struct
{
  size_t row;
  size_t column;
  double value;
} DblNonzero;

int compareDblNonzeros(const void* pa, const void* pb)
{
  size_t aRow = ((const DblNonzero*)pa)->row;
  size_t bRow = ((const DblNonzero*)pb)->row;
  if (aRow < bRow)
    return -1;
  else if (aRow > bRow)
    return +1;

  int aColumn = ((const DblNonzero*)pa)->column;
  int bColumn = ((const DblNonzero*)pb)->column;
  return aColumn - bColumn;
}

CMR_ERROR CMRdblmatCreateFromSparseStream(CMR* cmr, FILE* stream, CMR_DBLMAT** presult)
{
  assert(cmr);
  assert(presult);
  assert(!*presult);
  assert(stream);

  size_t numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%lu %lu %lu", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return CMR_ERROR_INPUT;

  /* Read all nonzeros. */

  DblNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
  size_t entry = 0;
  for (size_t i = 0; i < numNonzeros; ++i)
  {
    size_t row;
    size_t column;
    double value;
    numRead = fscanf(stream, "%lu %lu %lf", &row, &column, &value);
    if (numRead < 3 || row == 0 || column == 0 || row > numRows || column > numColumns)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      return CMR_ERROR_INPUT;
    }
    if (value != 0.0)
    {
      nonzeros[entry].row = row - 1;
      nonzeros[entry].column = column - 1;
      nonzeros[entry].value = value;
      ++entry;
    }
  }
  numNonzeros = entry;

  /* We sort all nonzeros by row and then by column. */
  CMR_CALL( CMRsort(cmr, numNonzeros, nonzeros, sizeof(DblNonzero), compareDblNonzeros) );

  CMR_CALL( CMRdblmatCreate(cmr, presult, numRows, numColumns, numNonzeros) );
  CMR_DBLMAT* result = *presult;
  size_t previousRow = SIZE_MAX;
  size_t previousColumn = SIZE_MAX;
  size_t* pentryColumn = result->entryColumns;
  double* pentryValue = result->entryValues;
  for (size_t entry = 0; entry < numNonzeros; ++entry)
  {
    size_t row = nonzeros[entry].row;
    size_t column = nonzeros[entry].column;
    if (row == previousRow && column == previousColumn)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      CMR_CALL( CMRdblmatFree(cmr, presult) );
      return CMR_ERROR_INPUT;
    }
    while (previousRow < row || previousRow == SIZE_MAX)
    {
      ++previousRow;
      result->rowSlice[previousRow] = entry;
    }
    *pentryColumn++ = nonzeros[entry].column;
    *pentryValue++ = nonzeros[entry].value;
    previousColumn = column;
  }
  while (previousRow < numRows || previousRow == SIZE_MAX)
  {
    ++previousRow;
    result->rowSlice[previousRow] = numNonzeros;
  }

  CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );

  return CMR_OKAY;
}

typedef struct
{
  size_t row;
  size_t column;
  int value;
} IntNonzero;

int compareIntNonzeros(const void* pa, const void* pb)
{
  size_t  aRow = ((const IntNonzero*)pa)->row;
  size_t  bRow = ((const IntNonzero*)pb)->row;
  if (aRow < bRow)
    return -1;
  else if (aRow > bRow)
    return +1;

  size_t  aColumn = ((const IntNonzero*)pa)->column;
  size_t  bColumn = ((const IntNonzero*)pb)->column;
  return aColumn - bColumn;
}

CMR_ERROR CMRintmatCreateFromSparseStream(CMR* cmr, FILE* stream, CMR_INTMAT** presult)
{
  assert(cmr);
  assert(presult);
  assert(!*presult);
  assert(stream);

  size_t numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%lu %lu %lu", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return CMR_ERROR_INPUT;

  /* Read all nonzeros. */

  IntNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
  size_t entry = 0;
  for (size_t i = 0; i < numNonzeros; ++i)
  {
    size_t row;
    size_t column;
    int value;
    numRead = fscanf(stream, "%lu %lu %d", &row, &column, &value);
    if (numRead < 3 || row == 0 || column == 0 || row > numRows || column > numColumns)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      return CMR_ERROR_INPUT;
    }
    if (value != 0)
    {
      nonzeros[entry].row = row - 1;
      nonzeros[entry].column = column - 1;
      nonzeros[entry].value = value;
      ++entry;
    }
  }
  numNonzeros = entry;

  /* We sort all nonzeros by row and then by column. */
  CMR_CALL( CMRsort(cmr, numNonzeros, nonzeros, sizeof(DblNonzero), compareIntNonzeros) );

  CMR_CALL( CMRintmatCreate(cmr, presult, numRows, numColumns, numNonzeros) );
  CMR_INTMAT* result = *presult;
  size_t previousRow = SIZE_MAX;
  size_t previousColumn = SIZE_MAX;
  size_t* pentryColumn = result->entryColumns;
  int* pentryValue = result->entryValues;
  for (size_t entry = 0; entry < numNonzeros; ++entry)
  {
    size_t row = nonzeros[entry].row;
    size_t column = nonzeros[entry].column;
    if (row == previousRow && column == previousColumn)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      CMR_CALL( CMRintmatFree(cmr, presult) );
      return CMR_ERROR_INPUT;
    }
    while (previousRow < row || previousRow == SIZE_MAX)
    {
      ++previousRow;
      result->rowSlice[previousRow] = entry;
    }
    *pentryColumn++ = nonzeros[entry].column;
    *pentryValue++ = nonzeros[entry].value;
    previousColumn = column;
  }
  while (previousRow < numRows || previousRow == SIZE_MAX)
  {
    ++previousRow;
    result->rowSlice[previousRow] = numNonzeros;
  }

  CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );

  return CMR_OKAY;
}

typedef struct
{
  size_t row;
  size_t column;
  char value;
} ChrNonzero;

int compareChrNonzeros(const void* pa, const void* pb)
{
  size_t  aRow = ((const ChrNonzero*)pa)->row;
  size_t  bRow = ((const ChrNonzero*)pb)->row;
  if (aRow < bRow)
    return -1;
  else if (aRow > bRow)
    return +1;

  size_t  aColumn = ((const ChrNonzero*)pa)->column;
  size_t  bColumn = ((const ChrNonzero*)pb)->column;
  return aColumn - bColumn;
}

CMR_ERROR CMRchrmatCreateFromSparseStream(CMR* cmr, FILE* stream, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(presult);
  assert(!*presult);
  assert(stream);

  size_t numRows, numColumns, numNonzeros;
  int numRead = fscanf(stream, "%lu %lu %lu", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
    return CMR_ERROR_INPUT;

  /* Read all nonzeros. */

  ChrNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
  size_t entry = 0;
  for (size_t i = 0; i < numNonzeros; ++i)
  {
    size_t row;
    size_t column;
    int value;
    numRead = fscanf(stream, "%lu %lu %d", &row, &column, &value);
    if (numRead < 3 || row == 0 || column == 0 || row > numRows || column > numColumns)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      return CMR_ERROR_INPUT;
    }
    if (value != 0)
    {
      nonzeros[entry].row = row - 1;
      nonzeros[entry].column = column - 1;
      nonzeros[entry].value = value;
      ++entry;
    }
  }
  numNonzeros = entry;

  /* We sort all nonzeros by row and then by column. */
  CMR_CALL( CMRsort(cmr, numNonzeros, nonzeros, sizeof(DblNonzero), compareIntNonzeros) );

  CMR_CALL( CMRchrmatCreate(cmr, presult, numRows, numColumns, numNonzeros) );
  CMR_CHRMAT* result = *presult;
  size_t previousRow = SIZE_MAX;
  size_t previousColumn = SIZE_MAX;
  size_t* pentryColumn = result->entryColumns;
  char* pentryValue = result->entryValues;
  for (size_t entry = 0; entry < numNonzeros; ++entry)
  {
    size_t row = nonzeros[entry].row;
    size_t column = nonzeros[entry].column;
    if (row == previousRow && column == previousColumn)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      CMR_CALL( CMRchrmatFree(cmr, presult) );
      return CMR_ERROR_INPUT;
    }
    while (previousRow < row || previousRow == SIZE_MAX)
    {
      ++previousRow;
      result->rowSlice[previousRow] = entry;
    }
    *pentryColumn++ = nonzeros[entry].column;
    *pentryValue++ = nonzeros[entry].value;
    previousColumn = column;
  }
  while (previousRow < numRows || previousRow == SIZE_MAX)
  {
    ++previousRow;
    result->rowSlice[previousRow] = numNonzeros;
  }

  CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatCreateFromDenseStream(CMR* cmr, FILE* stream, CMR_DBLMAT** presult)
{
  assert(cmr);
  assert(presult);
  assert(!*presult);
  assert(stream);

  size_t numRows, numColumns;
  int numRead = fscanf(stream, "%lu %lu", &numRows, &numColumns);
  if (numRead < 2)
    return CMR_ERROR_INPUT;

  CMR_CALL( CMRdblmatCreate(cmr, presult, numRows, numColumns, 0) );
  CMR_DBLMAT* result = *presult;

  /* Initial memory. */
  size_t memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  size_t* entryColumns = NULL;
  double* entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &entryColumns, memEntries) );
  CMR_CALL( CMRallocBlockArray(cmr, &entryValues, memEntries) );

  size_t entry = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    result->rowSlice[row] = entry;
    for (size_t column = 0; column < numColumns; ++column)
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
  result->rowSlice[numRows] = entry;

  /* Make arrays smaller again. */
  if (entry < memEntries)
  {
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry) );
  }

  result->entryColumns = entryColumns;
  result->entryValues = entryValues;
  result->numNonzeros = entry;

  return CMR_OKAY;
}

CMR_ERROR CMRintmatCreateFromDenseStream(CMR* cmr, FILE* stream, CMR_INTMAT** presult)
{
  assert(cmr);
  assert(presult);
  assert(!*presult);
  assert(stream);

  size_t numRows, numColumns;
  int numRead = fscanf(stream, "%lu %lu", &numRows, &numColumns);
  if (numRead < 2)
    return CMR_ERROR_INPUT;

  CMR_CALL( CMRintmatCreate(cmr, presult, numRows, numColumns, 0) );
  CMR_INTMAT* result = *presult;

  /* Initial memory. */
  size_t memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  size_t* entryColumns = NULL;
  int* entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &entryColumns, memEntries) );
  CMR_CALL( CMRallocBlockArray(cmr, &entryValues, memEntries) );

  size_t entry = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    result->rowSlice[row] = entry;
    for (size_t column = 0; column < numColumns; ++column)
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
  result->rowSlice[numRows] = entry;

  /* Make arrays smaller again. */
  if (entry < memEntries)
  {
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry) );
  }

  result->entryColumns = entryColumns;
  result->entryValues = entryValues;
  result->numNonzeros = entry;

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatCreateFromDenseStream(CMR* cmr, FILE* stream, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(presult);
  assert(!*presult);
  assert(stream);

  size_t numRows, numColumns;
  int numRead = fscanf(stream, "%lu %lu", &numRows, &numColumns);
  if (numRead < 2)
    return CMR_ERROR_INPUT;

  CMR_CALL( CMRchrmatCreate(cmr, presult, numRows, numColumns, 0) );
  CMR_CHRMAT* result = *presult;

  /* Initial memory. */
  size_t memEntries = numRows * numColumns;
  if (memEntries > 256)
    memEntries = 256;
  size_t* entryColumns = NULL;
  char* entryValues = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &entryColumns, memEntries) );
  CMR_CALL( CMRallocBlockArray(cmr, &entryValues, memEntries) );

  size_t entry = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    result->rowSlice[row] = entry;
    for (size_t column = 0; column < numColumns; ++column)
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
  result->rowSlice[numRows] = entry;

  /* Make arrays smaller again. */
  if (entry < memEntries)
  {
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry) );
  }

  result->entryColumns = entryColumns;
  result->entryValues = entryValues;
  result->numNonzeros = entry;

  return CMR_OKAY;
}

bool CMRdblmatCheckEqual(CMR_DBLMAT* matrix1, CMR_DBLMAT* matrix2)
{
  CMRconsistencyAssert( CMRdblmatConsistency(matrix1) );
  CMRconsistencyAssert( CMRdblmatConsistency(matrix2) );

  if (matrix1->numRows != matrix2->numRows)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;
  if (matrix1->numNonzeros != matrix2->numNonzeros)
    return false;

  for (size_t row = 0; row < matrix1->numRows; ++row)
  {
    size_t first1 = matrix1->rowSlice[row];
    size_t first2 = matrix2->rowSlice[row];
    if (first1 != first2)
      return false;
    size_t beyond1 = matrix1->rowSlice[row + 1];
    size_t beyond2 = matrix2->rowSlice[row + 1];
    if (beyond1 != beyond2)
      return false;

    for (size_t entry = first1; entry < beyond1; ++entry)
    {
      if (matrix1->entryColumns[entry] != matrix2->entryColumns[entry])
        return false;
      if (matrix1->entryValues[entry] != matrix2->entryValues[entry])
        return false;
    }
  }

  return true;
}

bool CMRintmatCheckEqual(CMR_INTMAT* matrix1, CMR_INTMAT* matrix2)
{
  CMRconsistencyAssert( CMRintmatConsistency(matrix1) );
  CMRconsistencyAssert( CMRintmatConsistency(matrix2) );

  if (matrix1->numRows != matrix2->numRows)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;
  if (matrix1->numNonzeros != matrix2->numNonzeros)
    return false;

  for (size_t row = 0; row < matrix1->numRows; ++row)
  {
    size_t first1 = matrix1->rowSlice[row];
    size_t first2 = matrix2->rowSlice[row];
    if (first1 != first2)
      return false;
    size_t beyond1 = matrix1->rowSlice[row + 1];
    size_t beyond2 = matrix2->rowSlice[row + 1];
    if (beyond1 != beyond2)
      return false;

    for (size_t entry = first1; entry < beyond1; ++entry)
    {
      if (matrix1->entryColumns[entry] != matrix2->entryColumns[entry])
        return false;
      if (matrix1->entryValues[entry] != matrix2->entryValues[entry])
        return false;
    }
  }

  return true;
}

bool CMRchrmatCheckEqual(CMR_CHRMAT* matrix1, CMR_CHRMAT* matrix2)
{
  CMRconsistencyAssert( CMRchrmatConsistency(matrix1) );
  CMRconsistencyAssert( CMRchrmatConsistency(matrix2) );

  if (matrix1->numRows != matrix2->numRows)
    return false;
  if (matrix1->numColumns != matrix2->numColumns)
    return false;
  if (matrix1->numNonzeros != matrix2->numNonzeros)
    return false;

  for (size_t row = 0; row < matrix1->numRows; ++row)
  {
    size_t first1 = matrix1->rowSlice[row];
    size_t first2 = matrix2->rowSlice[row];
    if (first1 != first2)
      return false;
    size_t beyond1 = matrix1->rowSlice[row + 1];
    size_t beyond2 = matrix2->rowSlice[row + 1];
    if (beyond1 != beyond2)
      return false;

    for (size_t entry = first1; entry < beyond1; ++entry)
    {
      if (matrix1->entryColumns[entry] != matrix2->entryColumns[entry])
        return false;
      if (matrix1->entryValues[entry] != matrix2->entryValues[entry])
        return false;
    }
  }

  return true;
}

CMR_ERROR CMRdblmatCheckTranspose(CMR* cmr, CMR_DBLMAT* matrix1, CMR_DBLMAT* matrix2, bool* pareTranspose)
{
  assert(cmr);
  CMRconsistencyAssert( CMRdblmatConsistency(matrix1) );
  CMRconsistencyAssert( CMRdblmatConsistency(matrix2) );
  assert(pareTranspose);

  if (matrix1->numRows != matrix2->numColumns || matrix1->numColumns != matrix2->numRows
    || matrix1->numNonzeros != matrix2->numNonzeros)
  {
    *pareTranspose = false;
    return CMR_OKAY;
  }

  *pareTranspose = true;
  size_t* currentColumn1 = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &currentColumn1, matrix2->numRows) );
  for (size_t row2 = 0; row2 < matrix2->numRows; ++row2)
    currentColumn1[row2] = matrix2->rowSlice[row2];

  for (size_t row1 = 0; row1 < matrix1->numRows; ++row1)
  {
    size_t first = matrix1->rowSlice[row1];
    size_t beyond = matrix1->rowSlice[row1 + 1];
    for (size_t entry1 = first; entry1 < beyond; ++entry1)
    {
      size_t column1 = matrix1->entryColumns[entry1];
      size_t entry2 = currentColumn1[column1];
      if (matrix2->entryColumns[entry2] != row1 || matrix2->entryValues[entry2] != matrix1->entryValues[entry1])
      {
        *pareTranspose = false;
        goto cleanup;
      }
      currentColumn1[column1]++;
    }
  }

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &currentColumn1) );

  free(currentColumn1);

  return CMR_OKAY;
}

CMR_ERROR CMRintmatCheckTranspose(CMR* cmr, CMR_INTMAT* matrix1, CMR_INTMAT* matrix2, bool* pareTranspose)
{
  assert(cmr);
  CMRconsistencyAssert( CMRintmatConsistency(matrix1) );
  CMRconsistencyAssert( CMRintmatConsistency(matrix2) );
  assert(pareTranspose);

  if (matrix1->numRows != matrix2->numColumns || matrix1->numColumns != matrix2->numRows
    || matrix1->numNonzeros != matrix2->numNonzeros)
  {
    *pareTranspose = false;
    return CMR_OKAY;
  }

  *pareTranspose = true;
  size_t* currentColumn1 = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &currentColumn1, matrix2->numRows) );
  for (size_t row2 = 0; row2 < matrix2->numRows; ++row2)
    currentColumn1[row2] = matrix2->rowSlice[row2];

  for (size_t row1 = 0; row1 < matrix1->numRows; ++row1)
  {
    size_t first = matrix1->rowSlice[row1];
    size_t beyond = matrix1->rowSlice[row1 + 1];
    for (size_t entry1 = first; entry1 < beyond; ++entry1)
    {
      size_t column1 = matrix1->entryColumns[entry1];
      size_t entry2 = currentColumn1[column1];
      if (matrix2->entryColumns[entry2] != row1 || matrix2->entryValues[entry2] != matrix1->entryValues[entry1])
      {
        *pareTranspose = false;
        goto cleanup;
      }
      currentColumn1[column1]++;
    }
  }

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &currentColumn1) );

  free(currentColumn1);

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatCheckTranspose(CMR* cmr, CMR_CHRMAT* matrix1, CMR_CHRMAT* matrix2, bool* pareTranspose)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix1) );
  CMRconsistencyAssert( CMRchrmatConsistency(matrix2) );
  assert(pareTranspose);

  if (matrix1->numRows != matrix2->numColumns || matrix1->numColumns != matrix2->numRows
    || matrix1->numNonzeros != matrix2->numNonzeros)
  {
    *pareTranspose = false;
    return CMR_OKAY;
  }

  *pareTranspose = true;
  size_t* currentColumn1 = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &currentColumn1, matrix2->numRows) );
  for (size_t row2 = 0; row2 < matrix2->numRows; ++row2)
    currentColumn1[row2] = matrix2->rowSlice[row2];

  for (size_t row1 = 0; row1 < matrix1->numRows; ++row1)
  {
    size_t first = matrix1->rowSlice[row1];
    size_t beyond = matrix1->rowSlice[row1 + 1];
    for (size_t entry1 = first; entry1 < beyond; ++entry1)
    {
      size_t column1 = matrix1->entryColumns[entry1];
      size_t entry2 = currentColumn1[column1];
      if (matrix2->entryColumns[entry2] != row1 || matrix2->entryValues[entry2] != matrix1->entryValues[entry1])
      {
        *pareTranspose = false;
        goto cleanup;
      }
      currentColumn1[column1]++;
    }
  }

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &currentColumn1) );

  free(currentColumn1);

  return CMR_OKAY;
}

char* CMRdblmatConsistency(CMR_DBLMAT* matrix)
{
  if (!matrix)
    return CMRconsistencyMessage("CMR_DBLMAT is NULL.");
  if (!matrix->rowSlice)
    return CMRconsistencyMessage("CMR_DBLMAT is does not have rowSlice array.");

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    int first = matrix->rowSlice[row];
    int beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      if (matrix->entryValues[entry] == 0.0)
      {
        return CMRconsistencyMessage("CMR_DBLMAT contains zero entry %lu in row %lu, column %lu.\n",
          entry, row, matrix->entryColumns[entry]);
      }
    }
    for (size_t entry = first + 1; entry < beyond; ++entry)
    {
      if (matrix->entryColumns[entry - 1] == matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_DBLMAT contains duplicate nonzeros in row %lu, column %lu, entries %lu and %lu.\n",
          row, matrix->entryColumns[entry], entry - 1, entry);
      }
      if (matrix->entryColumns[entry - 1] > matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_DBLMAT contains nonzeros in row %lu in wrong order: entry %lu has column %lu and entry %lu has column %lu.\n",
          row, entry - 1, matrix->entryColumns[entry - 1], entry, matrix->entryColumns[entry]);
      }
    }
  }

  return NULL;
}

char* CMRintmatConsistency(CMR_INTMAT* matrix)
{
  if (!matrix)
    return CMRconsistencyMessage("CMR_INTMAT is NULL.");
  if (!matrix->rowSlice)
    return CMRconsistencyMessage("CMR_INTMAT is does not have rowSlice array.");

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    int first = matrix->rowSlice[row];
    int beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      if (matrix->entryValues[entry] == 0)
      {
        return CMRconsistencyMessage("CMR_INTMAT contains zero entry %lu in row %lu, column %lu.\n",
          entry, row, matrix->entryColumns[entry]);
      }
    }
    for (size_t entry = first + 1; entry < beyond; ++entry)
    {
      if (matrix->entryColumns[entry - 1] == matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_INTMAT contains duplicate nonzeros in row %lu, column %lu, entries %lu and %lu.\n",
          row, matrix->entryColumns[entry], entry - 1, entry);
      }
      if (matrix->entryColumns[entry - 1] > matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_INTMAT contains nonzeros in row %lu in wrong order: entry %lu has column %lu and entry %lu has column %lu.\n",
          row, entry - 1, matrix->entryColumns[entry - 1], entry, matrix->entryColumns[entry]);
      }
    }
  }

  return NULL;
}

char* CMRchrmatConsistency(CMR_CHRMAT* matrix)
{
  if (!matrix)
    return CMRconsistencyMessage("CMR_CHRMAT is NULL.");
  if (!matrix->rowSlice)
    return CMRconsistencyMessage("CMR_CHRMAT is does not have rowSlice array.");

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    int first = matrix->rowSlice[row];
    int beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      if (matrix->entryValues[entry] == 0)
      {
        return CMRconsistencyMessage("CMR_CHRMAT contains zero entry %lu in row %lu, column %lu.\n",
          entry, row, matrix->entryColumns[entry]);
      }
    }
    for (size_t entry = first + 1; entry < beyond; ++entry)
    {
      if (matrix->entryColumns[entry - 1] == matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_CHRMAT contains duplicate nonzeros in row %lu, column %lu, entries %lu and %lu.\n",
          row, matrix->entryColumns[entry], entry - 1, entry);
      }
      if (matrix->entryColumns[entry - 1] > matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_CHRMAT contains nonzeros in row %lu in wrong order: entry %lu has column %lu and entry %lu has column %lu.\n",
          row, entry - 1, matrix->entryColumns[entry - 1], entry, matrix->entryColumns[entry]);
      }
    }
  }

  return NULL;
}

bool CMRdblmatIsBinary(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  CMRconsistencyAssert( CMRdblmatConsistency(matrix) );
  assert(psubmatrix || !*psubmatrix);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      double value = matrix->entryValues[entry];
      double rounded = round(value);
      if (rounded < 0.0 || rounded > 1.0 || fabs(value - rounded) > epsilon)
      {
        if (psubmatrix)
          CMR_CALL( CMRsubmatCreate1x1(cmr, row, matrix->entryColumns[entry], psubmatrix) );
        return false;
      }
    }
  }

  return true;
}

bool CMRintmatIsBinary(CMR* cmr, CMR_INTMAT* matrix, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );
  assert(psubmatrix || !*psubmatrix);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      int value = matrix->entryValues[entry];
      if (value < 0 || value > 1)
      {
        if (psubmatrix)
          CMR_CALL( CMRsubmatCreate1x1(cmr, row, matrix->entryColumns[entry], psubmatrix) );
        return false;
      }
    }
  }

  return true;
}

bool CMRchrmatIsBinary(CMR* cmr, CMR_CHRMAT* matrix, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(psubmatrix || !*psubmatrix);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      int value = matrix->entryValues[entry];
      if (value < 0 || value > 1)
      {
        if (psubmatrix)
          CMR_CALL( CMRsubmatCreate1x1(cmr, row, matrix->entryColumns[entry], psubmatrix) );
        return false;
      }
    }
  }

  return true;
}

bool CMRdblmatIsTernary(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  CMRconsistencyAssert( CMRdblmatConsistency(matrix) );
  assert(!psubmatrix || !*psubmatrix);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      double value = matrix->entryValues[entry];
      double rounded = round(value);
      if (rounded < -1.0 || rounded > 1.0 || fabs(value - rounded) > epsilon)
      {
        if (psubmatrix)
          CMR_CALL( CMRsubmatCreate1x1(cmr, row, matrix->entryColumns[entry], psubmatrix) );
        return false;
      }
    }
  }

  return true;
}

bool CMRintmatIsTernary(CMR* cmr, CMR_INTMAT* matrix, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );
  assert(!psubmatrix || !*psubmatrix);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      int value = matrix->entryValues[entry];
      if (value < -1 || value > 1)
      {
        if (psubmatrix)
          CMR_CALL( CMRsubmatCreate1x1(cmr, row, matrix->entryColumns[entry], psubmatrix) );
        return false;
      }
    }
  }

  return true;
}

bool CMRchrmatIsTernary(CMR* cmr, CMR_CHRMAT* matrix, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(!psubmatrix || !*psubmatrix);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      int value = matrix->entryValues[entry];
      if (value < -1 || value > 1)
      {
        if (psubmatrix)
          CMR_CALL( CMRsubmatCreate1x1(cmr, row, matrix->entryColumns[entry], psubmatrix) );
        return false;
      }
    }
  }

  return true;
}

CMR_ERROR CMRdblmatFilterSubmat(CMR* cmr, CMR_DBLMAT* matrix, CMR_SUBMAT* submatrix, CMR_DBLMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRdblmatConsistency(matrix) );
  assert(submatrix);
  assert(presult && !*presult);

  int* columnMap = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnMap, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnMap[column] = -1;
  for (size_t j = 0; j < submatrix->numColumns; ++j)
  {
    assert(submatrix->columns[j] < matrix->numColumns);
    columnMap[submatrix->columns[j]] = j;
  }

  /* Count nonzeros. */
  size_t numNonzeros = 0;
  for (size_t i = 0; i < submatrix->numRows; ++i)
  {
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    size_t first = matrix->rowSlice[r];
    size_t beyond = matrix->rowSlice[r + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t c = matrix->entryColumns[entry];
      if (columnMap[c] >= 0)
        ++numNonzeros;
    }
  }

  CMR_CALL( CMRdblmatCreate(cmr, presult, submatrix->numRows, submatrix->numColumns, numNonzeros) );
  CMR_DBLMAT* result = *presult;

  /* Copy nonzeros. */
  result->numNonzeros = 0;
  for (size_t i = 0; i < submatrix->numRows; ++i)
  {
    result->rowSlice[i] = result->numNonzeros;
    size_t r = submatrix->rows[i];
    assert(r < matrix->numRows);

    size_t first = matrix->rowSlice[r];
    size_t beyond = matrix->rowSlice[r + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t column = matrix->entryColumns[entry];
      if (columnMap[column] >= 0)
      {
        result->entryColumns[result->numNonzeros] = columnMap[column];
        result->entryValues[result->numNonzeros] = matrix->entryValues[entry];
        result->numNonzeros++;
      }
    }
  }
  result->rowSlice[result->numRows] = result->numNonzeros;
  CMR_CALL( CMRfreeStackArray(cmr, &columnMap) );

  /* Sort the rows. */
  CMR_CALL( CMRdblmatSortNonzeros(cmr, result) );

  CMRconsistencyAssert( CMRdblmatConsistency(result) );

  return CMR_OKAY;
}

CMR_ERROR CMRintmatFilterSubmat(CMR* cmr, CMR_INTMAT* matrix, CMR_SUBMAT* submatrix, CMR_INTMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );
  assert(submatrix);
  assert(presult && !*presult);

  int* columnMap = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnMap, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnMap[column] = -1;
  for (size_t j = 0; j < submatrix->numColumns; ++j)
  {
    assert(submatrix->columns[j] < matrix->numColumns);
    columnMap[submatrix->columns[j]] = j;
  }

  /* Count nonzeros. */
  size_t numNonzeros = 0;
  for (size_t i = 0; i < submatrix->numRows; ++i)
  {
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    size_t first = matrix->rowSlice[r];
    size_t beyond = matrix->rowSlice[r + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t c = matrix->entryColumns[entry];
      if (columnMap[c] >= 0)
        ++numNonzeros;
    }
  }

  CMR_CALL( CMRintmatCreate(cmr, presult, submatrix->numRows, submatrix->numColumns, numNonzeros) );
  CMR_INTMAT* result = *presult;

  /* Copy nonzeros. */
  result->numNonzeros = 0;
  for (size_t i = 0; i < submatrix->numRows; ++i)
  {
    result->rowSlice[i] = result->numNonzeros;
    size_t r = submatrix->rows[i];
    assert(r < matrix->numRows);

    size_t first = matrix->rowSlice[r];
    size_t beyond = matrix->rowSlice[r + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t column = matrix->entryColumns[entry];
      if (columnMap[column] >= 0)
      {
        result->entryColumns[result->numNonzeros] = columnMap[column];
        result->entryValues[result->numNonzeros] = matrix->entryValues[entry];
        result->numNonzeros++;
      }
    }
  }
  result->rowSlice[result->numRows] = result->numNonzeros;
  CMR_CALL( CMRfreeStackArray(cmr, &columnMap) );

  /* Sort the rows. */
  CMR_CALL( CMRintmatSortNonzeros(cmr, result) );

  CMRconsistencyAssert( CMRintmatConsistency(result) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatFilter(CMR* cmr, CMR_CHRMAT* matrix, size_t numRows, size_t* rows, size_t numColumns,
  size_t* columns, CMR_CHRMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(presult && !*presult);

  int* columnMap = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnMap, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnMap[column] = -1;
  for (size_t j = 0; j < numColumns; ++j)
  {
    assert(columns[j] < matrix->numColumns);
    columnMap[columns[j]] = j;
  }

  /* Count nonzeros. */
  size_t numNonzeros = 0;
  for (size_t i = 0; i < numRows; ++i)
  {
    int r = rows[i];
    assert(r < matrix->numRows);

    size_t first = matrix->rowSlice[r];
    size_t beyond = matrix->rowSlice[r + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t c = matrix->entryColumns[entry];
      if (columnMap[c] >= 0)
        ++numNonzeros;
    }
  }

  CMR_CALL( CMRchrmatCreate(cmr, presult, numRows, numColumns, numNonzeros) );
  CMR_CHRMAT* result = *presult;

  /* Copy nonzeros. */
  result->numNonzeros = 0;
  for (size_t i = 0; i < numRows; ++i)
  {
    result->rowSlice[i] = result->numNonzeros;
    size_t r = rows[i];
    assert(r < matrix->numRows);

    size_t first = matrix->rowSlice[r];
    size_t beyond = matrix->rowSlice[r + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t column = matrix->entryColumns[entry];
      if (columnMap[column] >= 0)
      {
        result->entryColumns[result->numNonzeros] = columnMap[column];
        result->entryValues[result->numNonzeros] = matrix->entryValues[entry];
        result->numNonzeros++;
      }
    }
  }
  result->rowSlice[result->numRows] = result->numNonzeros;
  CMR_CALL( CMRfreeStackArray(cmr, &columnMap) );

  /* Sort the rows. */
  CMR_CALL( CMRchrmatSortNonzeros(cmr, result) );

  CMRconsistencyAssert( CMRchrmatConsistency(result) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatFilterSubmat(CMR* cmr, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix, CMR_CHRMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(submatrix);
  assert(presult && !*presult);

  CMR_CALL( CMRchrmatFilter(cmr, matrix, submatrix->numRows, submatrix->rows, submatrix->numColumns,
    submatrix->columns, presult) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatSupport(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(!*presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  size_t resultEntry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    result->rowSlice[row] = resultEntry;
    int first = matrix->rowSlice[row];
    int beyond = matrix->rowSlice[row + 1];
    for (size_t matrixEntry = first; matrixEntry < beyond; ++matrixEntry)
    {
      if (fabs(matrix->entryValues[matrixEntry]) > epsilon)
      {
        result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
        result->entryValues[resultEntry] = 1;
        resultEntry++;
      }
    }
  }
  result->rowSlice[matrix->numRows] = resultEntry;
  CMR_CALL( CMRchrmatChangeNumNonzeros(cmr, result, resultEntry) );

  return CMR_OKAY;
}

CMR_ERROR CMRintmatSupport(CMR* cmr, CMR_INTMAT* matrix, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(!*presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  size_t resultEntry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    result->rowSlice[row] = resultEntry;
    int first = matrix->rowSlice[row];
    int beyond = matrix->rowSlice[row + 1];
    for (size_t matrixEntry = first; matrixEntry < beyond; ++matrixEntry)
    {
      result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
      result->entryValues[resultEntry] = 1;
      resultEntry++;
    }
  }
  result->rowSlice[matrix->numRows] = resultEntry;

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatSupport(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(!*presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  size_t resultEntry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    result->rowSlice[row] = resultEntry;
    int first = matrix->rowSlice[row];
    int beyond = matrix->rowSlice[row + 1];
    for (size_t matrixEntry = first; matrixEntry < beyond; ++matrixEntry)
    {
      result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
      result->entryValues[resultEntry] = 1;
      resultEntry++;
    }
  }
  result->rowSlice[matrix->numRows] = resultEntry;

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatSignedSupport(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, CMR_CHRMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRdblmatConsistency(matrix) );
  assert(epsilon >= 0);
  assert(presult);
  assert(!*presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  size_t resultEntry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    result->rowSlice[row] = resultEntry;
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t matrixEntry = first; matrixEntry < beyond; ++matrixEntry)
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
  result->rowSlice[matrix->numRows] = resultEntry;
  CMR_CALL( CMRchrmatChangeNumNonzeros(cmr, result, resultEntry) );

  return CMR_OKAY;
}

CMR_ERROR CMRintmatSignedSupport(CMR* cmr, CMR_INTMAT* matrix, CMR_CHRMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );
  assert(presult);
  assert(!*presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  size_t resultEntry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    result->rowSlice[row] = resultEntry;
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t matrixEntry = first; matrixEntry < beyond; ++matrixEntry)
    {
      result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
      result->entryValues[resultEntry] = (matrix->entryValues[matrixEntry] > 0) ? 1 : -1;
      resultEntry++;
    }
  }
  result->rowSlice[matrix->numRows] = resultEntry;

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatSignedSupport(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(presult);
  assert(!*presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  size_t resultEntry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    result->rowSlice[row] = resultEntry;
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t matrixEntry = first; matrixEntry < beyond; ++matrixEntry)
    {
      result->entryColumns[resultEntry] = matrix->entryColumns[matrixEntry];
      result->entryValues[resultEntry] = (matrix->entryValues[matrixEntry] > 0) ? 1 : -1;
      resultEntry++;
    }
  }
  result->rowSlice[matrix->numRows] = resultEntry;

  return CMR_OKAY;
}

CMR_ERROR CMRintmatToChr(CMR* cmr, CMR_INTMAT* matrix, CMR_CHRMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );
  assert(presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  for (size_t row = 0; row <= matrix->numRows; ++row)
    result->rowSlice[row] = matrix->rowSlice[row];

  for (size_t e = 0; e < matrix->numNonzeros; ++e)
  {
    result->entryColumns[e] = matrix->entryColumns[e];
    int x = matrix->entryValues[e];
    if (x <= CHAR_MAX && x >= CHAR_MIN)
      result->entryValues[e] = x;
    else
      return CMR_ERROR_OVERFLOW;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatFindEntry(CMR* cmr, CMR_DBLMAT* matrix, size_t row, size_t column, size_t* pentry)
{
  assert(cmr);
  CMRconsistencyAssert( CMRdblmatConsistency(matrix) );
  assert(pentry);

  size_t lower = matrix->rowSlice[row];
  size_t upper = matrix->rowSlice[row + 1];
  while (lower < upper)
  {
    size_t entry = (lower + upper) / 2;
    size_t c = matrix->entryColumns[entry];
    if (c < column)
      lower = entry + 1;
    else if (c > column)
      upper = entry;
    else
    {
      *pentry = entry;
      return CMR_OKAY;
    }
  }
  *pentry = SIZE_MAX;

  return CMR_OKAY;
}

CMR_ERROR CMRintmatFindEntry(CMR* cmr, CMR_INTMAT* matrix, size_t row, size_t column, size_t* pentry)
{
  assert(cmr);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );
  assert(pentry);

  size_t lower = matrix->rowSlice[row];
  size_t upper = matrix->rowSlice[row + 1];
  while (lower < upper)
  {
    size_t entry = (lower + upper) / 2;
    size_t c = matrix->entryColumns[entry];
    if (c < column)
      lower = entry + 1;
    else if (c > column)
      upper = entry;
    else
    {
      *pentry = entry;
      return CMR_OKAY;
    }
  }
  *pentry = SIZE_MAX;

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatFindEntry(CMR* cmr, CMR_CHRMAT* matrix, size_t row, size_t column, size_t* pentry)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(pentry);

  size_t lower = matrix->rowSlice[row];
  size_t upper = matrix->rowSlice[row + 1];
  while (lower < upper)
  {
    size_t entry = (lower + upper) / 2;
    size_t c = matrix->entryColumns[entry];
    if (c < column)
      lower = entry + 1;
    else if (c > column)
      upper = entry;
    else
    {
      *pentry = entry;
      return CMR_OKAY;
    }
  }
  *pentry = SIZE_MAX;

  return CMR_OKAY;
}

// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/matrix.h>

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "sort.h"
#include "env_internal.h"
#include "listmatrix.h"

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
  CMRfreeBlock(cmr, psubmatrix);

  return CMR_OKAY;
}

CMR_ERROR CMRsubmatTranspose(CMR_SUBMAT* submatrix)
{
  assert(submatrix);

  size_t tempSize = submatrix->numRows;
  submatrix->numRows = submatrix->numColumns;
  submatrix->numColumns = tempSize;

  size_t* tempArray = submatrix->rows;
  submatrix->rows = submatrix->columns;
  submatrix->columns = tempArray;

  return CMR_OKAY;
}

CMR_ERROR CMRsubmatZoomSubmat(CMR* cmr, CMR_SUBMAT* reference, CMR_SUBMAT* input, CMR_SUBMAT** poutput)
{
  assert(cmr);
  assert(reference);
  assert(input);
  assert(poutput);

  CMR_CALL( CMRsubmatCreate(cmr, input->numRows, input->numColumns, poutput) );
  CMR_SUBMAT* output = *poutput;

  /* Create reverse row mapping. */
  size_t numRows = 0;
  for (size_t r = 0; r < reference->numRows; ++r)
  {
    size_t row = reference->rows[r];
    numRows = row > numRows ? row : numRows;
  }
  ++numRows;
  size_t* reverseRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &reverseRows, numRows) );
  for (size_t row = 0; row < numRows; ++row)
    reverseRows[row] = SIZE_MAX;
  for (size_t r = 0; r < reference->numRows; ++r)
    reverseRows[reference->rows[r]] = r;

  /* Create reverse column mapping. */
  size_t numColumns = 0;
  for (size_t c = 0; c < reference->numColumns; ++c)
  {
    size_t column = reference->columns[c];
    numColumns = column > numColumns ? column : numColumns;
  }
  ++numColumns;
  size_t* reverseColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &reverseColumns, numColumns) );
  for (size_t column = 0; column < numColumns; ++column)
    reverseColumns[column] = SIZE_MAX;
  for (size_t c = 0; c < reference->numColumns; ++c)
    reverseColumns[reference->columns[c]] = c;

  /* Fill submatrix. */
  for (size_t r = 0; r < input->numRows; ++r)
  {
    size_t row = input->rows[r];
    size_t submatrixRow = reverseRows[row];
    if (submatrixRow == SIZE_MAX)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &reverseColumns) );
      CMR_CALL( CMRfreeStackArray(cmr, &reverseRows) );
      CMR_CALL( CMRsubmatFree(cmr, poutput) );
      return CMR_ERROR_INPUT;
    }
    output->rows[r] = submatrixRow;
  }
  for (size_t c = 0; c < input->numColumns; ++c)
  {
    size_t column = input->columns[c];
    size_t submatrixColumn = reverseColumns[column];
    if (submatrixColumn == SIZE_MAX)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &reverseColumns) );
      CMR_CALL( CMRfreeStackArray(cmr, &reverseRows) );
      CMR_CALL( CMRsubmatFree(cmr, poutput) );
      return CMR_ERROR_INPUT;
    }
    output->columns[c] = submatrixColumn;
  }

  /* Cleanup. */

  CMR_CALL( CMRfreeStackArray(cmr, &reverseColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &reverseRows) );

  return CMR_OKAY;
}

CMR_ERROR CMRsubmatPrint(CMR* cmr, CMR_SUBMAT* submatrix, size_t numRows, size_t numColumns, FILE* stream)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(submatrix);
  assert(stream);

  fprintf(stream, "%zu %zu %zu %zu\n", numRows, numColumns, submatrix->numRows, submatrix->numColumns);
  for (size_t row = 0; row < submatrix->numRows; ++row)
    fprintf(stream, "%zu ", submatrix->rows[row] + 1);
  fputc('\n', stream);
  for (size_t column = 0; column < submatrix->numColumns; ++column)
    fprintf(stream, "%zu ", submatrix->columns[column] + 1);
  fputc('\n', stream);

  return CMR_OKAY;
}

CMR_ERROR CMRsubmatWriteToFile(CMR* cmr, CMR_SUBMAT* submatrix, size_t numRows, size_t numColumns, const char* fileName)
{
  assert(cmr);
  assert(submatrix);

  FILE* stream;
  if (!fileName || !strcmp(fileName, "-"))
    stream = stdout;
  else
  {
    stream = fopen(fileName, "w");
    if (!stream)
      return CMR_ERROR_OUTPUT;
  }

  CMR_CALL( CMRsubmatPrint(cmr, submatrix, numRows, numColumns, stream) );

  if (stream != stdout)
    fclose(stream);

  return CMR_OKAY;
}

CMR_ERROR CMRsubmatReadFromStream(CMR* cmr, CMR_SUBMAT ** psubmatrix, size_t* pnumMatrixRows, size_t* pnumMatrixColumns,
  FILE* stream)
{
  assert(cmr);
  assert(psubmatrix);
  assert(stream);

  size_t numOriginalRows;
  size_t numOriginalColumns;
  size_t numRows;
  size_t numColumns;
  if (fscanf(stream, "%zu %zu %zu %zu", &numOriginalRows, &numOriginalColumns, &numRows, &numColumns) != 4)
    return CMR_ERROR_INPUT;

  if (numRows > numOriginalRows || numColumns > numOriginalColumns)
    return CMR_ERROR_INPUT;

  CMR_CALL( CMRsubmatCreate(cmr, numRows, numColumns, psubmatrix) );
  CMR_SUBMAT* submatrix = *psubmatrix;
  for (size_t r = 0; r < numRows; ++r)
  {
    size_t row;
    if (fscanf(stream, "%zu", &row) != 1)
      return CMR_ERROR_INPUT;

    if (row == 0 || row > numRows)
      return CMR_ERROR_INPUT;

    submatrix->rows[r] = row - 1;
  }
  for (size_t c = 0; c < numColumns; ++c)
  {
    size_t column;
    if (fscanf(stream, "%zu", &column) != 1)
      return CMR_ERROR_INPUT;

    if (column == 0 || column > numColumns)
      return CMR_ERROR_INPUT;

    submatrix->columns[c] = column - 1;
  }

  if (pnumMatrixRows)
    *pnumMatrixRows = numOriginalRows;
  if (pnumMatrixColumns)
    *pnumMatrixColumns = numOriginalColumns;

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
  if (matrix->entryColumns)
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
  if (matrix->entryColumns)
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
  if (matrix->entryColumns)
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
  assert(matrix);

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
    CMRdbgMsg(2, "Sorting nonzero entries in range [%zu,%zu).\n", first, beyond);
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

CMR_ERROR CMRdblmatPermute(CMR* cmr, CMR_DBLMAT* matrix, size_t* rows, size_t* columns, CMR_DBLMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);

  CMR_CALL( CMRdblmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_DBLMAT* result = *presult;

  size_t* columnsToResultColumns = NULL;
  if (columns)
  {
    CMR_CALL( CMRallocStackArray(cmr, &columnsToResultColumns, matrix->numColumns) );
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columnsToResultColumns[columns[column]] = column;
  }

  size_t resultEntry = 0;
  for (size_t resultRow = 0; resultRow < result->numRows; ++resultRow)
  {
    result->rowSlice[resultRow] = resultEntry;

    size_t row = rows ? rows[resultRow] : resultRow;
    CMRdbgMsg(0, "New row %zu is old row %zu.\n", resultRow, row);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t e = first; e < beyond; ++e)
    {
      result->entryValues[resultEntry] = matrix->entryValues[e];
      result->entryColumns[resultEntry] =
        columnsToResultColumns ? columnsToResultColumns[matrix->entryColumns[e]] : matrix->entryColumns[e];
      CMRdbgMsg(2, "Entry in old column %zu is now in column %zu; value = %d\n", matrix->entryColumns[e],
        result->entryColumns[resultEntry], result->entryValues[resultEntry]);
      ++resultEntry;
    }
  }
  result->rowSlice[result->numRows] = resultEntry;

  if (columnsToResultColumns) 
    CMR_CALL( CMRfreeStackArray(cmr, &columnsToResultColumns) );

  CMR_CALL( CMRdblmatSortNonzeros(cmr, result) );

  return CMR_OKAY;
}

CMR_ERROR CMRintmatPermute(CMR* cmr, CMR_INTMAT* matrix, size_t* rows, size_t* columns, CMR_INTMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);

  CMR_CALL( CMRintmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_INTMAT* result = *presult;

  size_t* columnsToResultColumns = NULL;
  if (columns)
  {
    CMR_CALL( CMRallocStackArray(cmr, &columnsToResultColumns, matrix->numColumns) );
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columnsToResultColumns[columns[column]] = column;
  }

  size_t resultEntry = 0;
  for (size_t resultRow = 0; resultRow < result->numRows; ++resultRow)
  {
    result->rowSlice[resultRow] = resultEntry;

    size_t row = rows ? rows[resultRow] : resultRow;
    CMRdbgMsg(0, "New row %zu is old row %zu.\n", resultRow, row);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t e = first; e < beyond; ++e)
    {
      result->entryValues[resultEntry] = matrix->entryValues[e];
      result->entryColumns[resultEntry] =
        columnsToResultColumns ? columnsToResultColumns[matrix->entryColumns[e]] : matrix->entryColumns[e];
      CMRdbgMsg(2, "Entry in old column %zu is now in column %zu; value = %d\n", matrix->entryColumns[e],
        result->entryColumns[resultEntry], result->entryValues[resultEntry]);
      ++resultEntry;
    }
  }
  result->rowSlice[result->numRows] = resultEntry;

  if (columnsToResultColumns) 
    CMR_CALL( CMRfreeStackArray(cmr, &columnsToResultColumns) );

  CMR_CALL( CMRintmatSortNonzeros(cmr, result) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatPermute(CMR* cmr, CMR_CHRMAT* matrix, size_t* rows, size_t* columns, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_CHRMAT* result = *presult;

  size_t* columnsToResultColumns = NULL;
  if (columns)
  {
    CMR_CALL( CMRallocStackArray(cmr, &columnsToResultColumns, matrix->numColumns) );
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columnsToResultColumns[columns[column]] = column;
  }

  size_t resultEntry = 0;
  for (size_t resultRow = 0; resultRow < result->numRows; ++resultRow)
  {
    result->rowSlice[resultRow] = resultEntry;

    size_t row = rows ? rows[resultRow] : resultRow;
    CMRdbgMsg(0, "New row %zu is old row %zu.\n", resultRow, row);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t e = first; e < beyond; ++e)
    {
      result->entryValues[resultEntry] = matrix->entryValues[e];
      result->entryColumns[resultEntry] =
        columnsToResultColumns ? columnsToResultColumns[matrix->entryColumns[e]] : matrix->entryColumns[e];
      CMRdbgMsg(2, "Entry in old column %zu is now in column %zu; value = %d\n", matrix->entryColumns[e],
        result->entryColumns[resultEntry], result->entryValues[resultEntry]);
      ++resultEntry;
    }
  }
  result->rowSlice[result->numRows] = resultEntry;

  if (columnsToResultColumns) 
    CMR_CALL( CMRfreeStackArray(cmr, &columnsToResultColumns) );

  CMR_CALL( CMRchrmatSortNonzeros(cmr, result) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatPrintSparse(CMR* cmr, CMR_DBLMAT* matrix, FILE* stream)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%zu %zu %zu\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t entry = first; entry < beyond; ++entry)
      fprintf(stream, "%zu %zu %g\n", row+1, matrix->entryColumns[entry]+1, matrix->entryValues[entry]);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRintmatPrintSparse(CMR* cmr, CMR_INTMAT* matrix, FILE* stream)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%zu %zu %zu\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t entry = first; entry < beyond; ++entry)
      fprintf(stream, "%zu %zu %d\n", row+1, matrix->entryColumns[entry]+1, matrix->entryValues[entry]);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatPrintSparse(CMR* cmr, CMR_CHRMAT* matrix, FILE* stream)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%zu %zu %zu\n\n", matrix->numRows, matrix->numColumns, matrix->numNonzeros);
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t entry = first; entry < beyond; ++entry)
      fprintf(stream, "%zu %zu %d\n", row+1, matrix->entryColumns[entry]+1, matrix->entryValues[entry]);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatPrintDense(CMR* cmr, CMR_DBLMAT* matrix, FILE* stream, char zeroChar, bool header)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%zu %zu\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%zu ", (column+1) % 10);
    fputs("\n  ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%zu| ", (row+1) % 10);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    size_t column = 0;
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t entryColumn = matrix->entryColumns[entry];
      for (; column < entryColumn; ++column)
        fprintf(stream, "%c ", zeroChar);
      fprintf(stream, "%g ", matrix->entryValues[entry]);
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
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%zu %zu\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%zu ", (column+1) % 10);
    fputs("\n  ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%zu| ", (row+1) % 10);
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
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(matrix);
  assert(stream);

  fprintf(stream, "%zu %zu\n", matrix->numRows, matrix->numColumns);
  if (header)
  {
    fputs("   ", stream);
#if defined(CMR_DEBUG)
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%2zu ", (column+1) % 10);
#else /* !CMR_DEBUG */
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fprintf(stream, "%zu ", (column+1) % 10);
#endif /* CMR_DEBUG */
    fputs("\n  ", stream);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%zu| ", (row+1) % 10);
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    size_t column = 0;
#if defined(CMR_DEBUG)
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t entryColumn = matrix->entryColumns[entry];
      for (; column < entryColumn; ++column)
        fprintf(stream, " %c ", zeroChar);
      fprintf(stream, "%2d ", matrix->entryValues[entry]);
      ++column;
    }
    for (; column < matrix->numColumns; ++column)
      fprintf(stream, " %c ", zeroChar);
#else /* !CMR_DEBUG */
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
#endif /* CMR_DEBUG */
    fputc('\n', stream);
  }
  fflush(stream);

  return CMR_OKAY;
}



typedef struct
{
  size_t row;
  size_t column;
  double value;
} DblNonzero;

static
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
  int numRead = fscanf(stream, "%zu %zu %zu", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
  {
    CMRraiseErrorMessage(cmr, "Could not read number of rows, columns and nonzeros.");
    return CMR_ERROR_INPUT;
  }

  /* Read all nonzeros. */

  DblNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
  size_t entry = 0;
  for (size_t i = 0; i < numNonzeros; ++i)
  {
    size_t row;
    size_t column;
    double value;
    numRead = fscanf(stream, "%zu %zu %lf", &row, &column, &value);
    if (numRead < 3 || row == 0 || column == 0 || row > numRows || column > numColumns)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      if (numRead == 2)
        CMRraiseErrorMessage(cmr, "Could not read a double value of nonzero #%zu.", entry);
      else
        CMRraiseErrorMessage(cmr, "Could not read nonzero #%zu.", entry);
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
      CMRraiseErrorMessage(cmr, "Duplicate nonzero at row %zu and column %zu.", row, column);
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

static
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
  int numRead = fscanf(stream, "%zu %zu %zu", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
  {
    CMRraiseErrorMessage(cmr, "Could not read number of rows, columns and nonzeros.");
    return CMR_ERROR_INPUT;
  }

  /* Read all nonzeros. */

  IntNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
  size_t entry = 0;
  for (size_t i = 0; i < numNonzeros; ++i)
  {
    size_t row;
    size_t column;
    int value;
    numRead = fscanf(stream, "%zu %zu %d", &row, &column, &value);
    if (numRead < 3 || row == 0 || column == 0 || row > numRows || column > numColumns)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      if (numRead == 2)
        CMRraiseErrorMessage(cmr, "Could not read an integer value of nonzero #%zu.", entry);
      else
        CMRraiseErrorMessage(cmr, "Could not read nonzero #%zu.", entry);
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
  CMR_CALL( CMRsort(cmr, numNonzeros, nonzeros, sizeof(IntNonzero), compareIntNonzeros) );

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
      CMRraiseErrorMessage(cmr, "Duplicate nonzero at row %zu and column %zu.", row, column);
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

static
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
  int numRead = fscanf(stream, "%zu %zu %zu", &numRows, &numColumns, &numNonzeros);
  if (numRead < 3)
  {
    CMRraiseErrorMessage(cmr, "Could not read number of rows, columns and nonzeros.");
    return CMR_ERROR_INPUT;
  }

  /* Read all nonzeros. */

  ChrNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nonzeros, numNonzeros) );
  size_t entry = 0;
  for (size_t i = 0; i < numNonzeros; ++i)
  {
    size_t row;
    size_t column;
    int value;
    numRead = fscanf(stream, "%zu %zu %d", &row, &column, &value);
    if (numRead < 3 || row == 0 || column == 0 || row > numRows || column > numColumns)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &nonzeros) );
      if (numRead == 2)
        CMRraiseErrorMessage(cmr, "Could not read an integer value of nonzero #%zu.", entry);
      else
        CMRraiseErrorMessage(cmr, "Could not read nonzero #%zu.", entry);
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
  CMR_CALL( CMRsort(cmr, numNonzeros, nonzeros, sizeof(ChrNonzero), compareChrNonzeros) );

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
      CMRraiseErrorMessage(cmr, "Duplicate nonzero at row %zu and column %zu.", row, column);
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

CMR_ERROR CMRdblmatCreateFromSparseFile(CMR* cmr, const char* fileName, const char* stdinName, CMR_DBLMAT** presult)
{
  FILE* inputFile = (!stdinName || strcmp(fileName, stdinName)) ? fopen(fileName, "r") : stdin;
  if (!inputFile)
  {
    CMRraiseErrorMessage(cmr, "Could not open file <%s>.", fileName);
    return CMR_ERROR_INPUT;
  }

  CMR_ERROR error = CMRdblmatCreateFromSparseStream(cmr, inputFile, presult);
  if (!error)
  {
    /* Attempt to read another token. */
    char token[16+4];
    size_t numRead = fscanf(inputFile, "%16s", token);
    if (numRead > 0 && strlen(token))
    {
      if (strlen(token) == 16)
        strcat(token, "...");
      CMRraiseErrorMessage(cmr, "Found unexpected token \"%s\" after having read a *sparse* %zux%zu matrix with %zu nonzeros.",
        token, (*presult)->numRows, (*presult)->numColumns, (*presult)->numNonzeros);
      CMRdblmatFree(cmr, presult);
      error = CMR_ERROR_INPUT;
    }
  }

  if (inputFile != stdin)
    fclose(inputFile);

  return error;
}

CMR_ERROR CMRintmatCreateFromSparseFile(CMR* cmr, const char* fileName, const char* stdinName, CMR_INTMAT** presult)
{
  FILE* inputFile = (!stdinName || strcmp(fileName, stdinName)) ? fopen(fileName, "r") : stdin;
  if (!inputFile)
  {
    CMRraiseErrorMessage(cmr, "Could not open file <%s>.", fileName);
    return CMR_ERROR_INPUT;
  }

  CMR_ERROR error = CMRintmatCreateFromSparseStream(cmr, inputFile, presult);
  if (!error)
  {
    /* Attempt to read another token. */
    char token[16+4];
    size_t numRead = fscanf(inputFile, "%16s", token);
    if (numRead > 0 && strlen(token))
    {
      if (strlen(token) == 16)
        strcat(token, "...");
      CMRraiseErrorMessage(cmr, "Found unexpected token \"%s\" after having read a *sparse* %zux%zu matrix with %zu nonzeros.",
        token, (*presult)->numRows, (*presult)->numColumns, (*presult)->numNonzeros);
      CMRintmatFree(cmr, presult);
      error = CMR_ERROR_INPUT;
    }
  }

  if (inputFile != stdin)
    fclose(inputFile);

  return error;
}

CMR_ERROR CMRchrmatCreateFromSparseFile(CMR* cmr, const char* fileName, const char* stdinName, CMR_CHRMAT** presult)
{
  FILE* inputFile = (!stdinName || strcmp(fileName, stdinName)) ? fopen(fileName, "r") : stdin;
  if (!inputFile)
  {
    CMRraiseErrorMessage(cmr, "Could not open file <%s>.", fileName);
    return CMR_ERROR_INPUT;
  }

  CMR_ERROR error = CMRchrmatCreateFromSparseStream(cmr, inputFile, presult);
  if (!error)
  {
    /* Attempt to read another token. */
    char token[16+4];
    size_t numRead = fscanf(inputFile, "%16s", token);
    if (numRead > 0 && strlen(token))
    {
      if (strlen(token) == 16)
        strcat(token, "...");
      CMRraiseErrorMessage(cmr, "Found unexpected token \"%s\" after having read a *sparse* %zux%zu matrix with %zu nonzeros.",
        token, (*presult)->numRows, (*presult)->numColumns, (*presult)->numNonzeros);
      CMRchrmatFree(cmr, presult);
      error = CMR_ERROR_INPUT;
    }
  }

  if (inputFile != stdin)
    fclose(inputFile);

  return error;
}

CMR_ERROR CMRdblmatCreateFromDenseStream(CMR* cmr, FILE* stream, CMR_DBLMAT** presult)
{
  assert(cmr);
  assert(presult);
  assert(!*presult);
  assert(stream);

  size_t numRows, numColumns;
  int numRead = fscanf(stream, "%zu %zu", &numRows, &numColumns);
  if (numRead < 2)
  {
    CMRraiseErrorMessage(cmr, "Could not read number of rows and columns.");
    return CMR_ERROR_INPUT;
  }

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
      {
        CMRraiseErrorMessage(cmr, "Could not read matrix entry in row %zu and column %zu.", row, column);
        return CMR_ERROR_INPUT;
      }

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
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry < 256 ? 256 : entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry < 256 ? 256 : entry) );
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
  int numRead = fscanf(stream, "%zu %zu", &numRows, &numColumns);
  if (numRead < 2)
  {
    CMRraiseErrorMessage(cmr, "Could not read number of rows and columns.");
    return CMR_ERROR_INPUT;
  }

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
      {
        CMRraiseErrorMessage(cmr, "Could not read matrix entry in row %zu and column %zu.", row, column);
        return CMR_ERROR_INPUT;
      }

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
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry < 256 ? 256 : entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry < 256 ? 256 : entry) );
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
  int numRead = fscanf(stream, "%zu %zu", &numRows, &numColumns);
  if (numRead < 2)
  {
    CMRraiseErrorMessage(cmr, "Could not read number of rows and columns.");
    return CMR_ERROR_INPUT;
  }

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
      {
        CMRraiseErrorMessage(cmr, "Could not read matrix entry in row %zu and column %zu.", row, column);
        return CMR_ERROR_INPUT;
      }

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
    CMR_CALL( CMRreallocBlockArray(cmr, &entryColumns, entry < 256 ? 256 : entry) );
    CMR_CALL( CMRreallocBlockArray(cmr, &entryValues, entry < 256 ? 256 : entry) );
  }

  result->entryColumns = entryColumns;
  result->entryValues = entryValues;
  result->numNonzeros = entry;

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatCreateFromDenseFile(CMR* cmr, const char* fileName, const char* stdinName, CMR_DBLMAT** presult)
{
  FILE* inputFile = (!stdinName || strcmp(fileName, stdinName)) ? fopen(fileName, "r") : stdin;
  if (!inputFile)
  {
    CMRraiseErrorMessage(cmr, "Could not open file <%s>.", fileName);
    return CMR_ERROR_INPUT;
  }

  CMR_ERROR error = CMRdblmatCreateFromDenseStream(cmr, inputFile, presult);
  if (!error)
  {
    /* Attempt to read another token. */
    char token[16+4];
    int numRead = fscanf(inputFile, "%16s", token);
    if (numRead != EOF && numRead > 0 && strlen(token))
    {
      if (strlen(token) == 16)
        strcat(token, "...");
      CMRraiseErrorMessage(cmr, "Found unexpected token \"%s\" after having read a *dense* %zux%zu matrix with %zu nonzeros.",
        token, (*presult)->numRows, (*presult)->numColumns, (*presult)->numNonzeros);
      CMRdblmatFree(cmr, presult);
      error = CMR_ERROR_INPUT;
    }
  }

  if (inputFile != stdin)
    fclose(inputFile);

  return error;
}

CMR_ERROR CMRintmatCreateFromDenseFile(CMR* cmr, const char* fileName, const char* stdinName, CMR_INTMAT** presult)
{
  FILE* inputFile = (!stdinName || strcmp(fileName, stdinName)) ? fopen(fileName, "r") : stdin;
  if (!inputFile)
  {
    CMRraiseErrorMessage(cmr, "Could not open file <%s>.", fileName);
    return CMR_ERROR_INPUT;
  }

  CMR_ERROR error = CMRintmatCreateFromDenseStream(cmr, inputFile, presult);
  if (!error)
  {
    /* Attempt to read another token. */
    char token[16+4];
    int numRead = fscanf(inputFile, "%16s", token);
    if (numRead != EOF && numRead > 0 && strlen(token))
    {
      if (strlen(token) == 16)
        strcat(token, "...");
      CMRraiseErrorMessage(cmr, "Found unexpected token \"%s\" after having read a *dense* %zux%zu matrix with %zu nonzeros.",
        token, (*presult)->numRows, (*presult)->numColumns, (*presult)->numNonzeros);
      CMRintmatFree(cmr, presult);
      error = CMR_ERROR_INPUT;
    }
  }

  if (inputFile != stdin)
    fclose(inputFile);

  return error;
}

CMR_ERROR CMRchrmatCreateFromDenseFile(CMR* cmr, const char* fileName, const char* stdinName, CMR_CHRMAT** presult)
{
  FILE* inputFile = (!stdinName || strcmp(fileName, stdinName)) ? fopen(fileName, "r") : stdin;
  if (!inputFile)
  {
    CMRraiseErrorMessage(cmr, "Could not open file <%s>.", fileName);
    return CMR_ERROR_INPUT;
  }

  CMR_ERROR error = CMRchrmatCreateFromDenseStream(cmr, inputFile, presult);
  if (!error)
  {
    /* Attempt to read another token. */
    char token[16+4];
    int numRead = fscanf(inputFile, "%16s", token);
    if (numRead != EOF && numRead > 0 && strlen(token))
    {
      if (strlen(token) == 16)
        strcat(token, "...");
      CMRraiseErrorMessage(cmr, "Found unexpected token \"%s\" after having read a *dense* %zux%zu matrix with %zu nonzeros.",
        token, (*presult)->numRows, (*presult)->numColumns, (*presult)->numNonzeros);
      CMRchrmatFree(cmr, presult);
      error = CMR_ERROR_INPUT;
    }
  }

  if (inputFile != stdin)
    fclose(inputFile);

  return error;
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

  return CMR_OKAY;
}

char* CMRdblmatConsistency(CMR_DBLMAT* matrix)
{
  if (!matrix)
    return CMRconsistencyMessage("CMR_DBLMAT is NULL.");
  if (!matrix->rowSlice)
    return CMRconsistencyMessage("CMR_DBLMAT is does not have rowSlice array.");
  if (matrix->rowSlice[matrix->numRows] != matrix->numNonzeros)
  {
    return CMRconsistencyMessage("CMR_DBLMAT has inconsistent last slice index (%zu) and #nonzeros (%zu)",
      matrix->rowSlice[matrix->numRows], matrix->numNonzeros);
  }

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first + 1; entry < beyond; ++entry)
    {
      if (matrix->entryColumns[entry - 1] == matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_DBLMAT contains duplicate nonzeros in row %zu, column %zu, entries %zu and %zu.\n",
          row, matrix->entryColumns[entry], entry - 1, entry);
      }
      if (matrix->entryColumns[entry - 1] > matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_DBLMAT contains nonzeros in row %zu in wrong order: entry #%zu has column %zu and entry #%zu has column %zu.\n",
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
  if (matrix->rowSlice[matrix->numRows] != matrix->numNonzeros)
  {
    return CMRconsistencyMessage("CMR_INTMAT has inconsistent last slice index (%zu) and #nonzeros (%zu)",
      matrix->rowSlice[matrix->numRows], matrix->numNonzeros);
  }

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      if (matrix->entryValues[entry] == 0)
      {
        return CMRconsistencyMessage("CMR_INTMAT contains zero entry #%zu in row %zu, column %zu.\n",
          entry, row, matrix->entryColumns[entry]);
      }
    }
    for (size_t entry = first + 1; entry < beyond; ++entry)
    {
      if (matrix->entryColumns[entry - 1] == matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_INTMAT contains duplicate nonzeros in row %zu, column %zu, entries %zu and %zu.\n",
          row, matrix->entryColumns[entry], entry - 1, entry);
      }
      if (matrix->entryColumns[entry - 1] > matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_INTMAT contains nonzeros in row %zu in wrong order: entry #%zu has column %zu and entry #%zu has column %zu.\n",
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
  if (matrix->rowSlice[matrix->numRows] != matrix->numNonzeros)
  {
    return CMRconsistencyMessage("CMR_CHRMAT has inconsistent last slice index (%zu) and #nonzeros (%zuf)",
      matrix->rowSlice[matrix->numRows], matrix->numNonzeros);
  }

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t entry = first; entry < beyond; ++entry)
    {
      if (matrix->entryValues[entry] == 0)
      {
        return CMRconsistencyMessage("CMR_CHRMAT contains zero entry #%zu in row %zu, column %zu.\n",
          entry, row, matrix->entryColumns[entry]);
      }
    }
    for (size_t entry = first + 1; entry < beyond; ++entry)
    {
      if (matrix->entryColumns[entry - 1] == matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_CHRMAT contains duplicate nonzeros in row %zu, column %zu, entries %zu and %zu.\n",
          row, matrix->entryColumns[entry], entry - 1, entry);
      }
      if (matrix->entryColumns[entry - 1] > matrix->entryColumns[entry])
      {
        return CMRconsistencyMessage("CMR_CHRMAT contains nonzeros in row %zu in wrong order: entry #%zu has column %zu and entry #%zu has column %zu.\n",
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
  assert(!psubmatrix || !*psubmatrix);

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

CMR_ERROR CMRdblmatZoomSubmat(CMR* cmr, CMR_DBLMAT* matrix, CMR_SUBMAT* submatrix, CMR_DBLMAT** presult)
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
    size_t r = submatrix->rows[i];
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

CMR_ERROR CMRintmatZoomSubmat(CMR* cmr, CMR_INTMAT* matrix, CMR_SUBMAT* submatrix, CMR_INTMAT** presult)
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
    size_t r = submatrix->rows[i];
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
    size_t r = rows[i];
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

CMR_ERROR CMRchrmatZoomSubmat(CMR* cmr, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix, CMR_CHRMAT** presult)
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
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
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
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
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
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
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

CMR_ERROR CMRchrmatToInt(CMR* cmr, CMR_CHRMAT* matrix, CMR_INTMAT** presult)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(presult);

  CMR_CALL( CMRintmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, matrix->numNonzeros) );
  CMR_INTMAT* result = *presult;

  for (size_t row = 0; row <= matrix->numRows; ++row)
    result->rowSlice[row] = matrix->rowSlice[row];

  for (size_t e = 0; e < matrix->numNonzeros; ++e)
  {
    result->entryColumns[e] = matrix->entryColumns[e];
    result->entryValues[e] = matrix->entryValues[e];
  }

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

CMR_ERROR CMRdblmatFindEntry(CMR_DBLMAT* matrix, size_t row, size_t column, size_t* pentry)
{
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

CMR_ERROR CMRintmatFindEntry(CMR_INTMAT* matrix, size_t row, size_t column, size_t* pentry)
{
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

CMR_ERROR CMRchrmatFindEntry(CMR_CHRMAT* matrix, size_t row, size_t column, size_t* pentry)
{
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

static
size_t findMaximum(size_t* array, size_t length, size_t* pmaxIndex)
{
  assert(length > 0);

  size_t maximum = array[0];
  *pmaxIndex = 0;

  for (size_t i = 1; i < length; ++i)
  {
    if (array[i] > maximum)
    {
      maximum = array[i];
      *pmaxIndex = i;
    }
  }
  return maximum;
}

static
CMR_ERROR findBadSubmatrixByMaximum(
  CMR* cmr,
  ListMat8* listmatrix,
  CMR_SUBMAT** psubmatrix
)
{
  size_t* rowNumBadEntries = NULL;
  size_t* columnNumBadEntries = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowNumBadEntries, listmatrix->numRows) );
  for (size_t row = 0; row < listmatrix->numRows; ++row)
    rowNumBadEntries[row] = 0;
  CMR_CALL( CMRallocStackArray(cmr, &columnNumBadEntries, listmatrix->numColumns) );
  for (size_t column = 0; column < listmatrix->numColumns; ++column)
    columnNumBadEntries[column] = 0;

  for (size_t row = 0; row < listmatrix->numRows; ++row)
  {
    for (ListMat8Nonzero* nonzero = listmatrix->rowElements[row].head.right;
      nonzero != &listmatrix->rowElements[row].head; nonzero = nonzero->right)
    {
      if (nonzero->special)
      {
        nonzero->special = 1;
        assert(row == nonzero->row);
        rowNumBadEntries[row]++;
        columnNumBadEntries[nonzero->column]++;
      }
    }
  }

  size_t rowMaximumIndex;
  size_t columnMaximumIndex;
  size_t numRemainingRows = listmatrix->numRows;
  size_t numRemainingColumns = listmatrix->numColumns;
  while (true)
  {
    size_t rowMaximum = findMaximum(rowNumBadEntries, listmatrix->numRows, &rowMaximumIndex);
    if (rowMaximum == 0)
      break;

    size_t columnMaximum = findMaximum(columnNumBadEntries, listmatrix->numColumns, &columnMaximumIndex);
    
    CMRdbgMsg(2, "row/column maxima are %zu and %zu\n", rowMaximum, columnMaximum);
    
    if (rowMaximum >= columnMaximum)
    {
      for (ListMat8Nonzero* nz = listmatrix->rowElements[rowMaximumIndex].head.right;
        nz != &listmatrix->rowElements[rowMaximumIndex].head; nz = nz->right)
      {
        CMRdbgMsg(4, "Removing nonzero at %zu,%zu with special %d.\n", nz->row, nz->column, nz->special);
        if (nz->special)
        {
          assert(columnNumBadEntries[nz->column] > 0);
          columnNumBadEntries[nz->column]--;
        }
        nz->above->below = nz->below;
        nz->below->above = nz->above;
      }
      rowNumBadEntries[rowMaximumIndex] = 0;
      listmatrix->rowElements[rowMaximumIndex].head.above->below = listmatrix->rowElements[rowMaximumIndex].head.below;
      listmatrix->rowElements[rowMaximumIndex].head.below->above = listmatrix->rowElements[rowMaximumIndex].head.above;
      numRemainingRows--;
    }
    else
    {
      for (ListMat8Nonzero* nz = listmatrix->columnElements[columnMaximumIndex].head.below;
        nz != &listmatrix->columnElements[columnMaximumIndex].head; nz = nz->below)
      {
        CMRdbgMsg(4, "Removing nonzero at %zu,%zu with special %d.\n", nz->row, nz->column, nz->special);
        if (nz->special)
        {
          assert(rowNumBadEntries[nz->row] > 0);
          rowNumBadEntries[nz->row]--;
        }
        nz->left->right = nz->right;
        nz->right->left = nz->left;
      }
      columnNumBadEntries[columnMaximumIndex] = 0;
      listmatrix->columnElements[columnMaximumIndex].head.left->right = listmatrix->columnElements[columnMaximumIndex].head.right;
      listmatrix->columnElements[columnMaximumIndex].head.right->left = listmatrix->columnElements[columnMaximumIndex].head.left;
      numRemainingColumns--;
    }
  }

  CMR_CALL( CMRsubmatCreate(cmr, numRemainingRows, numRemainingColumns, psubmatrix) );
  CMR_SUBMAT* submatrix = *psubmatrix;
  numRemainingRows = 0;
  numRemainingColumns = 0;
  for (ListMat8Nonzero* rowHead = listmatrix->anchor.below; rowHead != &listmatrix->anchor; rowHead = rowHead->below)
  {
    submatrix->rows[numRemainingRows++] = rowHead->row;
  }
  for (ListMat8Nonzero* columnHead = listmatrix->anchor.right; columnHead != &listmatrix->anchor;
    columnHead = columnHead->right)
  {
    submatrix->columns[numRemainingColumns++] = columnHead->column;
  }

  CMR_CALL( CMRfreeStackArray(cmr, &columnNumBadEntries) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowNumBadEntries) );
  CMR_CALL( CMRlistmat8Free(cmr, &listmatrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatFindBinarySubmatrix(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(matrix);
  assert(epsilon >= 0.0);
  assert(psubmatrix);

  ListMat8* listmatrix = NULL;
  CMR_CALL( CMRlistmat8Alloc(cmr, matrix->numRows, matrix->numColumns, matrix->numNonzeros, &listmatrix) );
  CMR_CALL( CMRlistmat8InitializeFromDoubleMatrix(cmr, listmatrix, matrix, epsilon) );

  CMRdbgMsg(2, "List matrix has %zu nonzeros.\n", listmatrix->numNonzeros);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    for (ListMat8Nonzero* nonzero = listmatrix->rowElements[row].head.right;
      nonzero != &listmatrix->rowElements[row].head; nonzero = nonzero->right)
    {
      if (nonzero->value < 0 || nonzero->value > 1)
        nonzero->special = 1;
    }
  }
  
  CMR_CALL( findBadSubmatrixByMaximum(cmr, listmatrix, psubmatrix) );

  return CMR_OKAY;
}

CMR_ERROR CMRdblmatFindTernarySubmatrix(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(matrix);
  assert(epsilon >= 0.0);
  assert(psubmatrix);

  ListMat8* listmatrix = NULL;
  CMR_CALL( CMRlistmat8Alloc(cmr, matrix->numRows, matrix->numColumns, matrix->numNonzeros, &listmatrix) );
  CMR_CALL( CMRlistmat8InitializeFromDoubleMatrix(cmr, listmatrix, matrix, epsilon) );

  CMRdbgMsg(2, "List matrix has %zu nonzeros.\n", listmatrix->numNonzeros);

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    for (ListMat8Nonzero* nonzero = listmatrix->rowElements[row].head.right;
      nonzero != &listmatrix->rowElements[row].head; nonzero = nonzero->right)
    {
      if (nonzero->value < -1 || nonzero->value > +1)
        nonzero->special = 1;
    }
  }
  
  CMR_CALL( findBadSubmatrixByMaximum(cmr, listmatrix, psubmatrix) );

  return CMR_OKAY;
}


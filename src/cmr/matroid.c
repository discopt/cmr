// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/matroid.h>

#include "env_internal.h"
#include "seymour_internal.h"
#include "matrix_internal.h"
#include "listmatrix.h"
#include "linear_algebra_internal.h"

#include <assert.h>
#include <string.h>

/**
 * \brief Carries out \p numPivots pivots on \p matrix and stores the result in \p *presult.
 *
 * Calculations are done modulo \p characteristic. If that value is negative then the pivots must be regular, i.e.,
 * over every field. If the result would differ over fields, \p *presult will be \c NULL and, if given, \p *pviolator
 * will refer to the submatrix indexed by an affected entry as well as indices of all previous pivots.
 */

static
CMR_ERROR computePivots(
  CMR* cmr,               /**< \ref CMR environment . */
  CMR_CHRMAT* matrix,     /**< Input matrix. */
  size_t numPivots,       /**< Number of pivots. */
  size_t* pivotRows,      /**< Row indices of pivots. */
  size_t* pivotColumns,   /**< Column indices of pivots. */
  int characteristic,     /**< Characteristic of field. */
  CMR_SUBMAT** pviolator, /**< Pointer for storing a violating (irregular) submatrix; may be \c NULL. */
  CMR_CHRMAT** presult    /**< Pointer for storing the resulting matrix. */
  )
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(!numPivots || pivotRows);
  assert(!numPivots || pivotColumns);

#if defined(CMR_DEBUG)
  CMRdbgMsg(2, "Applying %zu pivot(s) to a %zux%zu matrix.\n", numPivots, matrix->numRows, matrix->numColumns);
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
  fflush(stdout);
#endif /* CMR_DEBUG */

  ListMat8* listmat = NULL;
  CMR_ERROR error = CMR_OKAY;
  size_t memNonzeros = 2 * matrix->numNonzeros + 256;
  CMR_CALL( CMRlistmat8Alloc(cmr, matrix->numRows, matrix->numColumns, memNonzeros, &listmat) );
  CMR_CALL( CMRlistmat8InitializeFromChrMatrix(cmr, listmat, matrix) );

  CMRdbgMsg(2, "#nonzeros = %zu, mem=%zu\n", matrix->numNonzeros, memNonzeros);

  size_t* affectedRows = NULL;
  size_t numAffectedRows;
  CMR_CALL( CMRallocStackArray(cmr, &affectedRows, matrix->numRows) );
  int* denseColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseColumn, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    denseColumn[row] = INT_MIN;

  size_t* affectedColumns = NULL;
  size_t numAffectedColumns;
  CMR_CALL( CMRallocStackArray(cmr, &affectedColumns, matrix->numColumns) );
  int* denseRow = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseRow, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    denseRow[column] = INT_MIN;

  for (size_t pivot = 0; pivot < numPivots; ++pivot)
  {
    size_t pivotRow = pivotRows[pivot];
    size_t pivotColumn = pivotColumns[pivot];

    CMRdbgMsg(4, "Applying pivot at r%zu,c%zu\n", pivotRow+1, pivotColumn+1);

    /* Compute the pivot value and all columns that are affected. */
    numAffectedColumns = 0;
    ListMat8Nonzero* head = &listmat->rowElements[pivotRow].head;
    for (ListMat8Nonzero* nz = head->right; nz != head; nz = nz->right)
    {
      if (denseRow[nz->column] == INT_MIN)
      {
        CMRdbgMsg(6, "Column c%zu is affected.\n", nz->column+1);
        affectedColumns[numAffectedColumns++] = nz->column;
        denseRow[nz->column] = 0;
      }
      denseRow[nz->column] += nz->value;
    }

    if (characteristic < 0)
    {
      /* Check for non-ternary entries in pivot row. */

      for (ListMat8Nonzero* nz = head->right; nz != head; nz = nz->right)
      {
        int value = denseRow[nz->column];
        if (value == INT_MIN || value >= -1 || value <= +1)
          continue;

        *presult = NULL;
        CMR_CALL( CMRsubmatCreate(cmr, pivot + 1, pivot + 1, pviolator) );
        CMR_SUBMAT* violator = *pviolator;
        for (size_t p = 0; p < pivot; ++p)
        {
          violator->rows[p] = pivotRows[p];
          violator->columns[p] = pivotColumns[p];
        }
        violator->rows[pivot] = nz->row;
        violator->columns[pivot] = nz->column;

        goto cleanup;
      }
    }

    int pivotValue = denseRow[pivotColumn];
    CMRdbgMsg(6, "Pivot value is %d.\n", pivotValue);
    if (pivotValue == INT_MIN || moduloNonnegative(pivotValue, characteristic) == 0)
    {
      error = CMR_ERROR_INPUT;
      goto cleanup;
    }

    /* Compute all rows that are affected. */
    head = &listmat->columnElements[pivotColumn].head;
    numAffectedRows = 0;
    for (ListMat8Nonzero* nz = head->below; nz != head; nz = nz->below)
    {
      if (denseColumn[nz->row] == INT_MIN)
      {
        CMRdbgMsg(6, "Row r%zu is affected.\n", nz->row+1);
        affectedRows[numAffectedRows++] = nz->row;
        denseColumn[nz->row] = 0;
      }
      denseColumn[nz->row] += nz->value;
    }

    if (characteristic < 0)
    {
      /* Check for non-ternary entries in pivot column. */

      for (ListMat8Nonzero* nz = head->below; nz != head; nz = nz->below)
      {
        int value = denseColumn[nz->row];
        if (value == INT_MIN || value >= -1 || value <= +1)
          continue;

        *presult = NULL;
        CMR_CALL( CMRsubmatCreate(cmr, pivot + 1, pivot + 1, pviolator) );
        CMR_SUBMAT* violator = *pviolator;
        for (size_t p = 0; p < pivot; ++p)
        {
          violator->rows[p] = pivotRows[p];
          violator->columns[p] = pivotColumns[p];
        }
        violator->rows[pivot] = nz->row;
        violator->columns[pivot] = nz->column;

        goto cleanup;
      }
    }

    /* Apply the pivot to all other entries first. */
    for (size_t r = 0; r < numAffectedRows; ++r)
    {
      size_t row = affectedRows[r];
      if (row == pivotRow)
        continue;

      int rowValue = denseColumn[row];
      if (rowValue == 0)
        continue;

      for (size_t c = 0; c < numAffectedColumns; ++c)
      {
        size_t column = affectedColumns[c];
        if (column == pivotColumn)
          continue;

        int columnValue = denseRow[column];
        if (columnValue == 0)
          continue;

        CMRdbgMsg(8, "Adding -%d*%d*%d = %d to r%zu,r%zu; #list nonzeros = %zu, mem list nonzeros = %zu\n", pivotValue,
          rowValue, columnValue, -pivotValue * rowValue * columnValue, row+1, column+1, listmat->numNonzeros,
          listmat->memNonzeros);

        CMR_CALL( CMRlistmat8Insert(cmr, listmat, row, column, -pivotValue * rowValue * columnValue, 0, NULL) );
      }
    }

#if defined(CMR_DEBUG)
    CMRdbgMsg(4, "Pivot at r%zu,c%zu; after updating non-pivot rows/columns:\n", pivotRow+1, pivotColumn+1);
    CMRlistmat8PrintDense(cmr, listmat, stdout);
#endif /* CMR_DEBUG */

    /* Apply the pivot to the pivot row and clear the dense vector. */
    head = &listmat->rowElements[pivotRow].head;
    for (ListMat8Nonzero* nz = head->right; nz != head; nz = nz->right)
    {
      if (pivotValue == -1 && nz->column != pivotColumn)
        nz->value *= -1;
      denseRow[nz->column] = INT_MIN;
    }

    /* Apply the pivot to the pivot column and clear the dense vector. */
    head = &listmat->columnElements[pivotColumn].head;
    for (ListMat8Nonzero* nz = head->below; nz != head; nz = nz->below)
    {
      if (pivotValue == -1 && nz->row != pivotRow)
        nz->value *= -1;
      else if (nz->row == pivotRow)
        nz->value *= -1;
      denseColumn[nz->row] = INT_MIN;
    }

#if defined(CMR_DEBUG)
    CMRdbgMsg(4, "Pivot at r%zu,c%zu; after updating all rows/columns:\n", pivotRow+1, pivotColumn+1);
    CMRlistmat8PrintDense(cmr, listmat, stdout);
#endif /* CMR_DEBUG */

  }

  /* Extract all nonzeros and merge duplicates. */
  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, memNonzeros) );
  CMR_CHRMAT* result = *presult;
  size_t entry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    /* Compute all columns with nonzeros in the row. */
    numAffectedColumns = 0;
    ListMat8Nonzero* head = &listmat->rowElements[row].head;
    for (ListMat8Nonzero* nz = head->right; nz != head; nz = nz->right)
    {
      CMRdbgMsg(6, "Considering nonzero %d at r%zu,c%zu.\n", nz->value, nz->row+1, nz->column+1);
      if (denseRow[nz->column] == INT_MIN)
      {
        affectedColumns[numAffectedColumns++] = nz->column;
        denseRow[nz->column] = 0;
      }
      denseRow[nz->column] += nz->value;
    }

    if (entry + numAffectedColumns > result->numNonzeros)
    {
      size_t newNumNonzeros = 2 * result->numNonzeros + 256;
      CMR_CALL( CMRchrmatChangeNumNonzeros(cmr, result, newNumNonzeros) );
    }

    CMRdbgMsg(8, "#entries in row r%zu is <= %zu.\n", row+1, numAffectedColumns);

    result->rowSlice[row] = entry;
    for (size_t i = 0; i < numAffectedColumns; ++i)
    {
      size_t column = affectedColumns[i];
      int value = denseRow[column];

      CMRdbgMsg(6, "Aggregated value in r%zu,c%zu is %d.\n", row+1, column+1, value);

      if (characteristic < 0 && (value < -1 || value > +1))
      {
        /* Found non-ternary entries in regular case. */
        CMR_CALL( CMRchrmatFree(cmr, presult) );

        CMR_CALL( CMRsubmatCreate(cmr, numPivots + 1, numPivots + 1, pviolator) );
        CMR_SUBMAT* violator = *pviolator;
        for (size_t p = 0; p < numPivots; ++p)
        {
          violator->rows[p] = pivotRows[p];
          violator->columns[p] = pivotColumns[p];
        }
        violator->rows[numPivots] = row;
        violator->columns[numPivots] = column;

        goto cleanup;
      }
      else
        value = moduloTernary(value, characteristic);

      denseRow[column] = INT_MIN;
      if (value == 0)
        continue;

      result->entryColumns[entry] = column;
      result->entryValues[entry] = value;
      ++entry;
    }
  }
  result->rowSlice[matrix->numRows] = entry;
  result->numNonzeros = entry;

  CMR_CALL( CMRchrmatSortNonzeros(cmr, result) );

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &denseRow) );
  CMR_CALL( CMRfreeStackArray(cmr, &affectedColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &denseColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &affectedRows) );
  CMR_CALL( CMRlistmat8Free(cmr, &listmat) );

  return error;
}

CMR_ERROR CMRchrmatBinaryPivot(CMR* cmr, CMR_CHRMAT* matrix, size_t pivotRow, size_t pivotColumn, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, 1, &pivotRow, &pivotColumn, 2, NULL, presult) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatTernaryPivot(CMR* cmr, CMR_CHRMAT* matrix, size_t pivotRow, size_t pivotColumn, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, 1, &pivotRow, &pivotColumn, 3, NULL, presult) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatRegularPivot(CMR* cmr, CMR_CHRMAT* matrix, size_t pivotRow, size_t pivotColumn,
  CMR_SUBMAT** pviolator, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, 1, &pivotRow, &pivotColumn, -3, pviolator, presult) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatBinaryPivots(CMR* cmr, CMR_CHRMAT* matrix, size_t numPivots, size_t* pivotRows, size_t* pivotColumns,
  CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(pivotRows);
  assert(pivotColumns);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, numPivots, pivotRows, pivotColumns, 2, NULL, presult) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatTernaryPivots(CMR* cmr, CMR_CHRMAT* matrix, size_t numPivots, size_t* pivotRows,
  size_t* pivotColumns, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(pivotRows);
  assert(pivotColumns);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, numPivots, pivotRows, pivotColumns, 3, NULL, presult) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatRegularPivots(CMR* cmr, CMR_CHRMAT* matrix, size_t numPivots, size_t* pivotRows,
  size_t* pivotColumns, CMR_SUBMAT** pviolator, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(pivotRows);
  assert(pivotColumns);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, numPivots, pivotRows, pivotColumns, -3, pviolator, presult) );

  return CMR_OKAY;
}


CMR_ERROR CMRminorCreate(CMR* cmr, CMR_MINOR** pminor, size_t numPivots, CMR_SUBMAT* submatrix, CMR_MINOR_TYPE type)
{
  assert(cmr);
  assert(pminor);
  assert(!*pminor);

  CMR_CALL( CMRallocBlock(cmr, pminor) );
  CMR_MINOR* minor = *pminor;
  minor->numPivots = numPivots;
  minor->pivotRows = NULL;
  minor->pivotColumns = NULL;
  if (numPivots)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &minor->pivotRows, numPivots) );
    CMR_CALL( CMRallocBlockArray(cmr, &minor->pivotColumns, numPivots) );
  }
  minor->remainingSubmatrix = submatrix;
  minor->type = type;

  return CMR_OKAY;
}

CMR_ERROR CMRminorFree(CMR* cmr, CMR_MINOR** pminor)
{
  assert(cmr);
  assert(pminor);

  if (!*pminor)
    return CMR_OKAY;

  CMR_MINOR* minor = *pminor;
  CMR_CALL( CMRfreeBlockArray(cmr, &minor->pivotRows) );
  CMR_CALL( CMRfreeBlockArray(cmr, &minor->pivotColumns) );
  CMR_CALL( CMRsubmatFree(cmr, &minor->remainingSubmatrix) );
  CMR_CALL( CMRfreeBlock(cmr, pminor) );

  return CMR_OKAY;
}

CMR_MINOR_TYPE CMRminorType(CMR_MINOR* minor)
{
  assert(minor);

  return minor->type;
}

size_t CMRminorNumPivots(CMR_MINOR* minor)
{
  assert(minor);

  return minor->numPivots;
}

size_t* CMRminorPivotRows(CMR_MINOR* minor)
{
  assert(minor);

  return minor->pivotRows;
}

size_t* CMRminorPivotColumns(CMR_MINOR* minor)
{
  assert(minor);

  return minor->pivotColumns;
}

CMR_SUBMAT* CMRminorSubmatrix(CMR_MINOR* minor)
{
  assert(minor);

  return minor->remainingSubmatrix;
}

CMR_ERROR CMRminorPrint(CMR* cmr, CMR_MINOR* minor, size_t numRows, size_t numColumns, FILE* stream)
{
  assert(cmr);
  assert(minor);
  assert(stream);

  CMR_CALL( CMRsubmatPrint(cmr, minor->remainingSubmatrix, numRows, numColumns, stream) );

  return CMR_OKAY;
}


CMR_ERROR CMRminorWriteToFile(CMR* cmr, CMR_MINOR* minor, size_t numRows, size_t numColumns, const char* fileName)
{
  assert(cmr);
  assert(minor);

  FILE* stream;
  if (!fileName || !strcmp(fileName, "-"))
    stream = stdout;
  else
  {
    stream = fopen(fileName, "w");
    if (!stream)
      return CMR_ERROR_OUTPUT;
  }

  CMR_CALL( CMRminorPrint(cmr, minor, numRows, numColumns, stream) );

  if (stream != stdout)
    fclose(stream);

  return CMR_OKAY;
}

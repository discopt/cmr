#include <tu/matrix.h>

#include <assert.h>
#include <stdlib.h>
#include <math.h>

void TUclearSparseDouble(TU_SPARSE_DOUBLE* sparse)
{
  assert(sparse);

  sparse->numRows = 0;
  sparse->numColumns = 0;
  sparse->numNonzeros = 0;

  free(sparse->rowStarts);
  sparse->rowStarts = NULL;

  free(sparse->entryColumns);
  sparse->entryColumns = NULL;

  free(sparse->entryValues);
  sparse->entryValues = NULL;
}

void TUclearSparseInt(TU_SPARSE_INT* sparse)
{
  assert(sparse);

  sparse->numRows = 0;
  sparse->numColumns = 0;
  sparse->numNonzeros = 0;

  free(sparse->rowStarts);
  sparse->rowStarts = NULL;

  free(sparse->entryColumns);
  sparse->entryColumns = NULL;

  free(sparse->entryValues);
  sparse->entryValues = NULL;
}

void TUclearSparseChar(TU_SPARSE_CHAR* sparse)
{
  assert(sparse);

  sparse->numRows = 0;
  sparse->numColumns = 0;
  sparse->numNonzeros = 0;

  if (sparse->rowStarts)
    free(sparse->rowStarts);
  sparse->rowStarts = NULL;

  if (sparse->entryColumns)
    free(sparse->entryColumns);
  sparse->entryColumns = NULL;

  if (sparse->entryValues)
    free(sparse->entryValues);
  sparse->entryValues = NULL;
}

void TUprintSparseAsDenseDouble(FILE* stream, TU_SPARSE_DOUBLE* sparse, char zeroChar, bool header)
{
  assert(stream != NULL);
  assert(sparse != NULL);
  double* rowEntries = (double*) calloc(sparse->numColumns, sizeof(double));

  fprintf(stream, "%d %d\n", sparse->numRows, sparse->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (int column = 0; column < sparse->numColumns; ++column)
      fprintf(stream, "%d ", column % 10);
    fputs("\n  ", stream);
    for (int column = 0; column < sparse->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (int row = 0; row < sparse->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%d| ", row % 10);
    int start = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int i = start; i < end; ++i)
      rowEntries[sparse->entryColumns[i]] = sparse->entryValues[i];
    for (int column = 0; column < sparse->numColumns; ++column)
    {
      double x = rowEntries[column];
      if (x == 0.0)
        fprintf(stream, "%c ", zeroChar);
      else
        fprintf(stream, "%f ", x);
    }
    for (int i = start; i < end; ++i)
      rowEntries[sparse->entryColumns[i]] = 0.0;
    fputc('\n', stream);
  }

  free(rowEntries);
}

void TUprintSparseAsDenseInt(FILE* stream, TU_SPARSE_INT* sparse, char zeroChar, bool header)
{
  assert(stream != NULL);
  assert(sparse != NULL);
  int* rowEntries = (int*) calloc(sparse->numColumns, sizeof(int));

  fprintf(stream, "%d %d\n", sparse->numRows, sparse->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (int column = 0; column < sparse->numColumns; ++column)
      fprintf(stream, "%d ", column % 10);
    fputs("\n  ", stream);
    for (int column = 0; column < sparse->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (int row = 0; row < sparse->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%d| ", row % 10);
    int start = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int i = start; i < end; ++i)
      rowEntries[sparse->entryColumns[i]] = sparse->entryValues[i];
    for (int column = 0; column < sparse->numColumns; ++column)
    {
      int x = rowEntries[column];
      if (x == 0.0)
        fprintf(stream, "%c ", zeroChar);
      else
        fprintf(stream, "%d ", x);
    }
    for (int i = start; i < end; ++i)
      rowEntries[sparse->entryColumns[i]] = 0.0;
    fputc('\n', stream);
  }

  free(rowEntries);
}

void TUprintSparseAsDenseChar(FILE* stream, TU_SPARSE_CHAR* sparse, char zeroChar, bool header)
{
  assert(stream != NULL);
  assert(sparse != NULL);
  char* rowEntries = (char*) calloc(sparse->numColumns, sizeof(char));

  fprintf(stream, "%d %d\n", sparse->numRows, sparse->numColumns);
  if (header)
  {
    fputs("   ", stream);
    for (int column = 0; column < sparse->numColumns; ++column)
      fprintf(stream, "%d ", column % 10);
    fputs("\n  ", stream);
    for (int column = 0; column < sparse->numColumns; ++column)
      fputs("--", stream);
    fputc('\n', stream);
  }
  for (int row = 0; row < sparse->numRows; ++row)
  {
    if (header)
      fprintf(stream, "%d| ", row % 10);
    int start = sparse->rowStarts[row];
    int end = row + 1 < sparse->numRows ? sparse->rowStarts[row + 1] : sparse->numNonzeros;
    for (int i = start; i < end; ++i)
      rowEntries[sparse->entryColumns[i]] = sparse->entryValues[i];
    for (int column = 0; column < sparse->numColumns; ++column)
    {
      char x = rowEntries[column];
      if (x == 0.0)
        fprintf(stream, "%c ", zeroChar);
      else
        fprintf(stream, "%d ", x);
    }
    for (int i = start; i < end; ++i)
      rowEntries[sparse->entryColumns[i]] = 0.0;
    fputc('\n', stream);
  }

  free(rowEntries);
}

bool TUcheckSparseEqualDouble(TU_SPARSE_DOUBLE* matrix1, TU_SPARSE_DOUBLE* matrix2)
{
  assert(TUcheckSparseSortedDouble(matrix1));
  assert(TUcheckSparseSortedDouble(matrix2));

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

bool TUcheckSparseEqualInt(TU_SPARSE_INT* matrix1, TU_SPARSE_INT* matrix2)
{
  assert(TUcheckSparseSortedInt(matrix1));
  assert(TUcheckSparseSortedInt(matrix2));

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

bool TUcheckSparseEqualChar(TU_SPARSE_CHAR* matrix1, TU_SPARSE_CHAR* matrix2)
{
  assert(TUcheckSparseSortedChar(matrix1));
  assert(TUcheckSparseSortedChar(matrix2));

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

bool TUcheckSparseTransposeDouble(TU_SPARSE_DOUBLE* matrix1, TU_SPARSE_DOUBLE* matrix2)
{
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
    for (int entry = begin; entry < end; ++entry)
    {
      int column = matrix1->entryColumns[entry];
      int entry2 = currentColumnEntries[column];
      if (matrix2->entryColumns[entry2] != row)
        return false;
      if (matrix2->entryValues[entry2] != matrix1->entryValues[entry])
        return false;
    }
  }

  free(currentColumnEntries);

  return true;
}

bool TUcheckSparseTransposeInt(TU_SPARSE_INT* matrix1, TU_SPARSE_INT* matrix2)
{
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
    for (int entry = begin; entry < end; ++entry)
    {
      int column = matrix1->entryColumns[entry];
      int entry2 = currentColumnEntries[column];
      if (matrix2->entryColumns[entry2] != row)
        return false;
      if (matrix2->entryValues[entry2] != matrix1->entryValues[entry])
        return false;
    }
  }

  free(currentColumnEntries);

  return true;
}

bool TUcheckSparseTransposeChar(TU_SPARSE_CHAR* matrix1, TU_SPARSE_CHAR* matrix2)
{
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
      if (matrix2->entryColumns[entry2] != row)
        return false;
      if (matrix2->entryValues[entry2] != matrix1->entryValues[entry1])
        return false;
      currentColumnEntries[column]++;
    }
  }

  free(currentColumnEntries);

  return true;
}

bool TUcheckSparseSortedDouble(TU_SPARSE_DOUBLE* sparse)
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

bool TUcheckSparseSortedInt(TU_SPARSE_INT* sparse)
{
  return TUcheckSparseSortedDouble((TU_SPARSE_DOUBLE*) sparse);
}

bool TUcheckSparseSortedChar(TU_SPARSE_CHAR* sparse)
{
  return TUcheckSparseSortedDouble((TU_SPARSE_DOUBLE*) sparse);
}

bool TUisTernaryDouble(TU_SPARSE_DOUBLE* sparse, double epsilon, TU_SUBMATRIX** submatrix)
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
          TUcreateSubmatrix1x1(submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool TUisTernaryInt(TU_SPARSE_INT* sparse, TU_SUBMATRIX** submatrix)
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
          TUcreateSubmatrix1x1(submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}

bool TUisTernaryChar(TU_SPARSE_CHAR* sparse, TU_SUBMATRIX** submatrix)
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
          TUcreateSubmatrix1x1(submatrix, row, sparse->entryColumns[entry]);
        return false;
      }
    }
  }

  return true;
}


void TUcreateSubmatrix(
  TU_SUBMATRIX** submatrix,
  int numRows,
  int numColumns
  )
{
  assert(submatrix != NULL);
  *submatrix = (TU_SUBMATRIX*) malloc(sizeof(TU_SUBMATRIX));
  (*submatrix)->numRows = numRows;
  (*submatrix)->rows = (int*) malloc(numRows * sizeof(int));
  (*submatrix)->numColumns = numColumns;
  (*submatrix)->columns = (int*) malloc(numColumns * sizeof(int));
}

void TUcreateSubmatrix1x1(
  TU_SUBMATRIX** submatrix, /**< Pointer to submatrix */
  int row, /**< Row of entry */
  int column /**< Column of entry */
  )
{
  TUcreateSubmatrix(submatrix, 1, 1);
  (*submatrix)->rows[0] = row;
  (*submatrix)->columns[0] = column;
}

void TUfreeSubmatrix(TU_SUBMATRIX** submatrix)
{
  assert(submatrix);

  if ((*submatrix)->rows)
    free((*submatrix)->rows);
  if ((*submatrix)->columns)
    free((*submatrix)->columns);
  free(*submatrix);
  *submatrix = NULL;
}

static int TUsortSubmatrixCompare(const void* p1, const void* p2)
{
  return *(int*)p1 - *(int*)p2;
}

void TUsortSubmatrix(TU_SUBMATRIX* submatrix)
{
  assert(submatrix);

  qsort(submatrix->rows, submatrix->numRows, sizeof(int), TUsortSubmatrixCompare);
  qsort(submatrix->columns, submatrix->numColumns, sizeof(int), TUsortSubmatrixCompare);
}

void TUfilterSubmatrixChar(TU_SPARSE_CHAR* matrix, TU_SUBMATRIX* submatrix, TU_SPARSE_CHAR* result)
{
  assert(matrix);
  assert(submatrix);
  assert(result);

  int* columnMap = (int*) malloc(matrix->numColumns * sizeof(int));
  for (int c = 0; c < matrix->numColumns; ++c)
    columnMap[c] = -1;
  for (int j = 0; j < submatrix->numColumns; ++j)
  {
    assert(submatrix->columns[j] < matrix->numColumns);
    columnMap[submatrix->columns[j]] = j;
  }

  result->numRows = submatrix->numRows;
  result->numColumns = submatrix->numColumns;
  result->numNonzeros = 0;
  result->rowStarts = (int*) malloc((submatrix->numRows + 1) * sizeof(int));

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
  
  result->entryColumns = (int*) malloc( numNonzeros * sizeof(int) );
  result->entryValues = (char*) malloc( numNonzeros * sizeof(char) );

  /* Copy nonzeros. */
  for (int i = 0; i < submatrix->numRows; ++i)
  {
    result->rowStarts[i] = result->numNonzeros;
    int r = submatrix->rows[i];
    assert(r < matrix->numRows);

    int begin = matrix->rowStarts[r];
    int end = r + 1 < matrix->numRows ? matrix->rowStarts[r+1] : matrix->numNonzeros;
    for (int e = begin; e < end; ++e)
    {
      int c = matrix->entryColumns[e];
      if (columnMap[c] >= 0)
      {
        result->entryColumns[result->numNonzeros] = columnMap[c];
        result->entryValues[result->numNonzeros] = matrix->entryValues[e];
        result->numNonzeros++;
      }
    }
  }
  result->rowStarts[result->numRows] = result->numNonzeros;

  if (columnMap)
    free(columnMap);
}

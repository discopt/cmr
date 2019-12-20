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

  free(sparse->rowStarts);
  sparse->rowStarts = NULL;

  free(sparse->entryColumns);
  sparse->entryColumns = NULL;

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

bool TUisTernaryDouble(TU_SPARSE_CHAR* sparse, double epsilon)
{
  assert(sparse != NULL);

  for (int entry = 0; entry < sparse->numNonzeros; ++entry)
  {
    double value = sparse->entryValues[entry];
    int rounded = (int)(value + 0.5);
    if (rounded < -1 || rounded > +1 || fabs(value - rounded) > epsilon)
      return false;
  }

  return true;
}

bool TUisTernaryInt(TU_SPARSE_INT* sparse)
{
  assert(sparse != NULL);

  for (int entry = 0; entry < sparse->numNonzeros; ++entry)
  {
    int value = sparse->entryValues[entry];
    if (value < -1 || value > +1)
      return false;
  }

  return true;
}

bool TUisTernaryChar(TU_SPARSE_CHAR* sparse)
{
  assert(sparse != NULL);

  for (int entry = 0; entry < sparse->numNonzeros; ++entry)
  {
    char value = sparse->entryValues[entry];
    if (value < -1 || value > +1)
      return false;
  }

  return true;
}


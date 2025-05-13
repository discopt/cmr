// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/named.h>

#include <stdint.h>
#include <float.h>

#include <cmr/camion.h>

static
void computeDegreeStats(
  CMR_CHRMAT* matrix,       /**< Input matrix. */
  size_t* rowDegrees,       /**< Array for storing degrees for row nodes. */
  size_t* columnDegrees,    /**< Array for storing degrees for column nodes. */
  size_t* rowDegreeStats,   /**< Array for storing how often each degree arises for row nodes. */
  size_t* columnDegreeStats /**< Array for storing how often each degree arises for column nodes. */
)
{
  /* Initialize column degrees to zero. */
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnDegrees[column] = 0;

  /* Compute row and column degrees. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = row + 1 < matrix->numRows ? matrix->rowSlice[row + 1] : matrix->numNonzeros;
    rowDegrees[row] = beyond - first;
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      columnDegrees[column]++;
    }
  }

  for (size_t rowDegree = 0; rowDegree <= matrix->numColumns; ++rowDegree)
    rowDegreeStats[rowDegree] = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
    rowDegreeStats[rowDegrees[row]]++;

  for (size_t columnDegree = 0; columnDegree <= matrix->numRows; ++columnDegree)
    columnDegreeStats[columnDegree] = 0;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnDegreeStats[columnDegrees[column]]++;
}

static
bool hasParallelSupports(
  CMR_CHRMAT* matrix,         /**< Input matrix. */
  size_t* rowSupportHash,     /**< Array for storing a hash per row. */
  size_t* columnSupportHash   /**< Array for storing a hash per column. */
)
{
  /* Initialize hashes as 0. */
  if (rowSupportHash)
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
      rowSupportHash[row] = 0;
  }
  if (columnSupportHash)
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columnSupportHash[column] = 0;
  }

  /* Update hashes. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = (row + 1 < matrix->numRows) ? matrix->rowSlice[row + 1] : matrix->numNonzeros;
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = matrix->entryColumns[e];
      if (rowSupportHash)
        rowSupportHash[row] += (1 << column);
      if (columnSupportHash)
        columnSupportHash[column] += (1 << row);
    }
  }

  /* Compare hashes pair-wise. */
  if (rowSupportHash)
  {
    for (size_t i = 0; i < matrix->numRows; ++i)
    {
      for (size_t j = i + 1; j < matrix->numRows; ++j)
      {
        if (rowSupportHash[i] == rowSupportHash[j])
          return true;
      }
    }
  }
  if (columnSupportHash)
  {
    for (size_t i = 0; i < matrix->numColumns; ++i)
    {
      for (size_t j = i + 1; j < matrix->numColumns; ++j)
      {
        if (columnSupportHash[i] == columnSupportHash[j])
          return true;
      }
    }
  }

  return false;
}

CMR_ERROR CMRisIdentityMatrix(CMR* cmr, CMR_CHRMAT* matrix, size_t* porder)
{
  assert(cmr);
  assert(matrix);
  assert(porder);
  CMR_UNUSED(cmr);

  *porder = SIZE_MAX;
  if (matrix->numRows != matrix->numColumns)
    return CMR_OKAY;

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = (row + 1 == matrix->numRows) ? matrix->numNonzeros : matrix->rowSlice[row+1];
    if ((beyond != first + 1) || (matrix->entryColumns[first] != first) || (matrix->entryValues[first] != 1))
      return CMR_OKAY;
  }

  *porder = matrix->numRows;
  return CMR_OKAY;
}

CMR_ERROR CMRcreateIdentityMatrix(CMR* cmr, size_t order, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(presult);

  CMR_CALL( CMRchrmatCreate(cmr, presult, order, order, order) );
  CMR_CHRMAT* matrix = *presult;
  for (size_t i = 0; i < order; ++i)
  {
    matrix->rowSlice[i] = i;
    matrix->entryColumns[i] = i;
    matrix->entryValues[i] = 1;
  }
  matrix->rowSlice[order] = order;

  return CMR_OKAY;
}

CMR_ERROR CMRisR10Matrix(CMR* cmr, CMR_CHRMAT* matrix, size_t* pisR10)
{
  assert(cmr);
  assert(matrix);
  assert(pisR10);

  if (matrix->numRows != 5 || matrix->numColumns != 5)
  {
    *pisR10 = 0;
    return CMR_OKAY;
  }

  size_t rowDegrees[5];
  size_t columnDegrees[5];
  size_t rowDegreeStats[6];
  size_t columnDegreeStats[6];
  computeDegreeStats(matrix, rowDegrees, columnDegrees, rowDegreeStats, columnDegreeStats);

  bool representative1 = (rowDegreeStats[3] == 4) && (rowDegreeStats[5] == 1) && (columnDegreeStats[3] == 4)
    && (columnDegreeStats[5] == 1);
  bool representative2 = (rowDegreeStats[3] == 5) && (columnDegreeStats[3] == 5);

  if (!representative1 && !representative2)
  {
    *pisR10 = 0;
    return CMR_OKAY;
  }

  size_t rowSupportHash[5];
  size_t columnSupportHash[5];
  bool parallel = hasParallelSupports(matrix, rowSupportHash, columnSupportHash);
  if (parallel)
  {
    *pisR10 = 0;
    return CMR_OKAY;
  }

  bool isSigned;
  CMR_CALL( CMRcamionTestSigns(cmr, matrix, &isSigned, NULL, NULL, DBL_MAX) );
  if (!isSigned)
    *pisR10 = 0;
  else
    *pisR10 = representative1 ? 1 : 2;

  return CMR_OKAY;
}


CMR_ERROR CMRcreateR10Matrix(CMR* cmr, size_t index, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(index >= 1 && index <= 2);
  assert(presult);

  if (index == 1)
  {
    size_t rowSlice[5] = { 0, 3, 6, 9, 12 };
    size_t entryColumns[17] = { 0, 3, 4,  0, 1, 4,  1, 2, 4,  2, 3, 4,  0, 1, 2, 3, 4 };

    CMR_CALL( CMRchrmatCreate(cmr, presult, 5, 5, 17) );
    CMR_CHRMAT* matrix = *presult;
    for (size_t row = 0; row < matrix->numRows; ++row)
      matrix->rowSlice[row] = rowSlice[row];
    for (size_t e = 0; e < matrix->numNonzeros; ++e)
    {
      matrix->entryColumns[e] = entryColumns[e];
      matrix->entryValues[e] = 1;
    }
  }
  else
  {
    size_t rowSlice[5] = { 0, 3, 6, 9, 12 };
    size_t entryColumns[15] = {  0, 1, 4,  0,  1, 2,  1,  2, 3,  2,  3, 4,  0, 3,  4 };
    char entryValues[15] =    { -1, 1, 1,  1, -1, 1,  1, -1, 1,  1, -1, 1,  1, 1, -1 };

    CMR_CALL( CMRchrmatCreate(cmr, presult, 5, 5, 15) );
    CMR_CHRMAT* matrix = *presult;
    for (size_t row = 0; row < matrix->numRows; ++row)
      matrix->rowSlice[row] = rowSlice[row];
    for (size_t e = 0; e < matrix->numNonzeros; ++e)
    {
      matrix->entryColumns[e] = entryColumns[e];
      matrix->entryValues[e] = entryValues[e];
    }
  }

  return CMR_OKAY;
}


CMR_ERROR CMRisR12Matrix(CMR* cmr, CMR_CHRMAT* matrix, size_t* pisR12)
{
  assert(cmr);
  assert(matrix);
  assert(pisR12);

  if (matrix->numRows != 6 || matrix->numColumns != 6)
  {
    *pisR12 = 0;
    return CMR_OKAY;
  }

  size_t rowDegrees[5];
  size_t columnDegrees[5];
  size_t rowDegreeStats[6];
  size_t columnDegreeStats[6];
  computeDegreeStats(matrix, rowDegrees, columnDegrees, rowDegreeStats, columnDegreeStats);



  size_t rowSupportHash[5];
  size_t columnSupportHash[5];
  bool parallel = hasParallelSupports(matrix, rowSupportHash, columnSupportHash);
  if (parallel)
  {
    *pisR12 = 0;
    return CMR_OKAY;
  }

  bool isSigned;
  CMR_CALL( CMRcamionTestSigns(cmr, matrix, &isSigned, NULL, NULL, DBL_MAX) );
  if (!isSigned)
    *pisR12 = 0;
  else
  {
    assert(false);
    return CMR_ERROR_INVALID;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRcreateR12Matrix(CMR* cmr, size_t index, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(index >= 1 && index <= 1);
  assert(presult);

  if (index == 1)
  {
    size_t rowSlice[6] = { 0, 3, 6, 10, 14, 17 };
    size_t entryColumns[20] = { 0, 2, 3,  1, 2, 3,  0, 2, 4, 5,  1,   3, 4, 5,  0, 2, 4,   1,  3, 5 };
    char entryValues[20] =    { 1, 1, 1,  1, 1, 1,  1, 1, 1, 1, -1,  -1, 1, 1,  1, 1, 1,  -1, -1, 1 };

    CMR_CALL( CMRchrmatCreate(cmr, presult, 6, 6, 20) );
    CMR_CHRMAT* matrix = *presult;
    for (size_t row = 0; row < matrix->numRows; ++row)
      matrix->rowSlice[row] = rowSlice[row];
    for (size_t e = 0; e < matrix->numNonzeros; ++e)
    {
      matrix->entryColumns[e] = entryColumns[e];
      matrix->entryValues[e] = entryValues[e];
    }
  }

  return CMR_OKAY;
}

CMR_ERROR CMRcreateK5Matrix(CMR* cmr, size_t index, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(index >= 1 && index <= 1);
  assert(presult);

  if (index == 1)
  {
    size_t rowSlice[4] = { 0, 3, 6, 9 };
    size_t entryColumns[12] = { 0, 3, 4,  0,  1,  5,  1, 2, 4,   2, 3, 5 };
    char entryValues[12] =    { 1, 1, 1,  1, -1, -1,  1, 1, 1,  -1, 1, 1 };

    CMR_CALL( CMRchrmatCreate(cmr, presult, 4, 6, 12) );
    CMR_CHRMAT* matrix = *presult;
    for (size_t row = 0; row < matrix->numRows; ++row)
      matrix->rowSlice[row] = rowSlice[row];
    for (size_t e = 0; e < matrix->numNonzeros; ++e)
    {
      matrix->entryColumns[e] = entryColumns[e];
      matrix->entryValues[e] = entryValues[e];
    }
  }

  return CMR_OKAY;
}

CMR_ERROR CMRcreateK33Matrix(CMR* cmr, size_t index, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(index >= 1 && index <= 1);
  assert(presult);

  if (index == 1)
  {
    size_t rowSlice[5] = { 0, 2, 4, 6, 8 };
    size_t entryColumns[12] = { 0, 1,  1, 2,  2, 3,  0, 3,  0, 1, 2, 3 };
    char entryValues[12] =    { 1, 1,  1, 1,  1, 1,  1, 1,  1, 1, 1, 1 };

    CMR_CALL( CMRchrmatCreate(cmr, presult, 5, 4, 12) );
    CMR_CHRMAT* matrix = *presult;
    for (size_t row = 0; row < matrix->numRows; ++row)
      matrix->rowSlice[row] = rowSlice[row];
    for (size_t e = 0; e < matrix->numNonzeros; ++e)
    {
      matrix->entryColumns[e] = entryColumns[e];
      matrix->entryValues[e] = entryValues[e];
    }
  }

  return CMR_OKAY;
}

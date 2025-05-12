// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/named.h>

#include <stdint.h>
#include <float.h>

#include <cmr/camion.h>


CMR_ERROR CMRisIdentityMatrix(CMR* cmr, CMR_CHRMAT* matrix, size_t* porder)
{
  assert(cmr);
  assert(matrix);
  assert(porder);

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

CMR_ERROR CMRisR10Matrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisR10)
{
  assert(cmr);
  assert(matrix);
  assert(pisR10);

  if (matrix->numRows != 5 || matrix->numColumns != 5)
    goto cleanup;

  size_t count3 = 0;
  size_t count5 = 0;
  for (size_t row = 0; row < 5; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = (row < 4) ? matrix->rowSlice[row + 1] : matrix->numNonzeros;
    size_t numNonzeros = beyond - first;
    if (numNonzeros == 3)
      ++count3;
    else if (numNonzeros == 5)
      ++count5;
    else
      goto cleanup;
  }

  if ((count3 != 5) && (count3 != 4 || count5 != 1))
    goto cleanup;

  /* No parallel rows. */
  size_t rowSupportHash[5] = { 0, 0, 0, 0, 0 };
  for (size_t row = 0; row < 5; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = (row < 4) ? matrix->rowSlice[row + 1] : matrix->numNonzeros;
    for (size_t e = first; e < beyond; ++e)
      rowSupportHash[row] += (1 << matrix->entryColumns[e]);
  }
  for (size_t i = 0; i < 5; ++i)
  {
    for (size_t j = i + 1; j < 5; ++j)
    {
      if (rowSupportHash[i] == rowSupportHash[j])
        goto cleanup;
    }
  }

  /* Column-wise nonzero counts. */
  size_t columnCounts[5] = { 0, 0, 0, 0, 0 };
  for (size_t e = 0; e < matrix->numNonzeros; ++e)
  {
    columnCounts[matrix->entryColumns[e]]++;
    if (matrix->entryValues[e] != 1 && matrix->entryValues[e] != -1)
      goto cleanup;
  }

  if (count3 == 5)
  {
    for (size_t c = 0; c < 5; ++c)
    {
      if (columnCounts[c] != 3)
        goto cleanup;
    }
  }
  else
  {
    count3 = 0;
    count5 = 0;
    for (size_t c = 0; c < 5; ++c)
    {
      if (columnCounts[c] == 3)
        count3++;
      else if (columnCounts[c] == 5)
        count5++;
    }

    if (count3 != 4 || count5 != 1)
      goto cleanup;

  }

  CMR_CALL( CMRcamionTestSigns(cmr, matrix, pisR10, NULL, NULL, DBL_MAX) );

  return CMR_OKAY;

cleanup:

  *pisR10 = false;
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
    size_t entryColumns[15] = {  0, 1, 4,  0,  1, 2,  1,  2, 3,  2,  3, 4,  3,  4, 0 };
    char entryValues[15] =    { -1, 1, 1,  1, -1, 1,  1, -1, 1,  1, -1, 1,  1, -1, 1 };

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

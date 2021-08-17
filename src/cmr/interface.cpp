#include "interface.h"

#include <cmr/env.h>

#include <cmr/total_unimodularity.hpp>

#include <boost/numeric/ublas/io.hpp>

extern "C"
CMR_ERROR CMRinterfaceTU(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTU, CMR_TU_DEC** pdec, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(matrix);
  assert(pisTU);
  assert(!psubmatrix || !*psubmatrix);

  tu::integer_matrix mat(matrix->numRows, matrix->numColumns, 0);
  for (size_t row = 0; row < (size_t)matrix->numRows; ++row)
  {
    size_t first = matrix->rowStarts[row];
    size_t beyond = (row + 1 < (size_t)matrix->numRows) ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (size_t i = first; i < beyond; ++i)
    {
      size_t column = matrix->entryColumns[i];
      mat(row,column) = matrix->entryValues[i];
    }
  }

  tu::submatrix_indices violator;
  if (psubmatrix)
  {
    *pisTU = tu::is_totally_unimodular(mat, violator);
  }
  else
  {
    *pisTU = tu::is_totally_unimodular(mat);
  }

  if (*pisTU && pdec)
    fprintf(stderr, "Retrieval of decomposition is not implemented, yet.");

  if (!violator.rows.empty())
  {
    CMR_CALL( CMRsubmatCreate(cmr, violator.rows.size(), violator.columns.size(), psubmatrix) );
    CMR_SUBMAT* submatrix = *psubmatrix;
    for (size_t row = 0; row < submatrix->numRows; ++row)
      submatrix->rows[row] = violator.rows[row];
    for (size_t column = 0; column < submatrix->numColumns; ++column)
      submatrix->columns[column] = violator.columns[column];
  }

  return CMR_OKAY;
}

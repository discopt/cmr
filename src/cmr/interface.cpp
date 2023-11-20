#include "interface.h"

#include <cmr/env.h>

#include "total_unimodularity.hpp"
#include "unimodularity.hpp"

#include <boost/numeric/ublas/io.hpp>

CMR_ERROR CMRinterfaceKModular(CMR* cmr, CMR_CHRMAT* matrix, size_t* pk)
{
  assert(cmr);
  assert(matrix);
  assert(pk);

  tu::integer_matrix mat(matrix->numRows, matrix->numColumns, 0);
  for (size_t row = 0; row < (size_t)matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t i = first; i < beyond; ++i)
    {
      size_t column = matrix->entryColumns[i];
      mat(row,column) = matrix->entryValues[i];
    }
  }

  size_t rank;
  unsigned int k;
  bool result = tu::is_k_modular(mat, rank, k);
  *pk = result ? k : 0;

  return CMR_OKAY;
}

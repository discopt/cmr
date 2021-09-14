#include "densematrix.h"

#include <assert.h>
#include <stdint.h>

#include "env_internal.h"

CMR_ERROR CMRdensebinmatrixCreateStack(CMR* cmr, size_t numRows, size_t numColumns, DenseBinaryMatrix** presult)
{
  assert(cmr);

  CMR_CALL( CMRallocStack(cmr, presult) );
  DenseBinaryMatrix* matrix = *presult;
  
  size_t size = (numRows * numColumns + (8 * sizeof(unsigned long long)) - 1) / (8 * sizeof(unsigned long long));
  matrix->numRows = numRows;
  matrix->numColumns = numColumns;
  matrix->data = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &matrix->data, size) );
  for (size_t i = 0; i < size; ++i)
    matrix->data[i] = 0UL;

  return CMR_OKAY;
}

CMR_ERROR CMRdensebinmatrixFreeStack(CMR* cmr, DenseBinaryMatrix** pmatrix)
{
  assert(cmr);

  CMR_CALL( CMRfreeStackArray(cmr, &(*pmatrix)->data) );
  CMR_CALL( CMRfreeStackArray(cmr, pmatrix) );

  return CMR_OKAY;
}

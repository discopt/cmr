// #define CMR_DEBUG /* Uncomment to debug this file. */

#include "env_internal.h"
#include "densematrix.h"

#include <assert.h>
#include <stdint.h>


CMR_ERROR CMRdensebinmatrixCreate(CMR* cmr, size_t numRows, size_t numColumns, DenseBinaryMatrix** presult)
{
  assert(cmr);

  CMR_CALL( CMRallocBlock(cmr, presult) );
  DenseBinaryMatrix* matrix = *presult;
  
  size_t size = (numRows * numColumns + (8 * sizeof(size_t)) - 1) / (8 * sizeof(size_t));
  CMRdbgMsg(10, "Creating %zux%zu DenseBinaryMatrix using %zu blocks each of which having %zu bytes.\n", numRows,
    numColumns, size, sizeof(size_t));
  matrix->numRows = numRows;
  matrix->numColumns = numColumns;
  matrix->data = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &matrix->data, size) );
  for (size_t i = 0; i < size; ++i)
    matrix->data[i] = (size_t) 0;

  return CMR_OKAY;
}

CMR_ERROR CMRdensebinmatrixFree(CMR* cmr, DenseBinaryMatrix** pmatrix)
{
  assert(cmr);
  assert(pmatrix);

  if (!*pmatrix)
    return CMR_OKAY;

  CMR_CALL( CMRfreeBlockArray(cmr, &(*pmatrix)->data) );
  CMR_CALL( CMRfreeBlockArray(cmr, pmatrix) );

  return CMR_OKAY;
}

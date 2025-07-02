#define CMR_DEBUG /* Uncomment to debug this file. */

#include "env_internal.h"
#include "densematrix.h"

#include <assert.h>
#include <stdint.h>


CMR_ERROR CMRdensebinmatrixCreate(CMR* cmr, size_t numRows, size_t numColumns, DenseBinaryMatrix** presult)
{
  assert(cmr);

  CMRdbgMsg(0, "sizeof(unsigned int) = %d\n", sizeof(unsigned int));
  CMRdbgMsg(0, "sizeof(unsigned long) = %d\n", sizeof(unsigned long));
  CMRdbgMsg(0, "sizeof(unsigned long int) = %d\n", sizeof(unsigned long int));
  CMRdbgMsg(0, "sizeof(unsigned long long) = %d\n", sizeof(unsigned long long));
  CMRdbgMsg(0, "sizeof(ssize_t) = %d\n", sizeof(ssize_t));
  CMRdbgMsg(0, "sizeof(size_t) = %d\n", sizeof(size_t));
  CMRdbgMsg(0, "sizeof(uint64_t) = %d\n", sizeof(uint64_t));
  CMRdbgMsg(0, "sizeof(uint32_t) = %d\n", sizeof(uint32_t));

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
    matrix->data[i] = 0UL;

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

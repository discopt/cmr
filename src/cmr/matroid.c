#include <cmr/matroid.h>

#include "env_internal.h"

#include <assert.h>
#include <string.h>

CMR_ERROR CMRminorCreate(CMR* cmr, CMR_MINOR** pminor, size_t numPivots, CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(pminor);
  assert(!*pminor);

  CMR_CALL( CMRallocBlock(cmr, pminor) );
  CMR_MINOR* minor = *pminor;
  minor->numPivots = numPivots;
  CMR_CALL( CMRallocBlockArray(cmr, &minor->pivotRows, numPivots) );
  CMR_CALL( CMRallocBlockArray(cmr, &minor->pivotColumns, numPivots) );
  minor->remainingSubmatrix = submatrix;

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

CMR_ERROR CMRminorWriteToStream(CMR* cmr, CMR_MINOR* minor, size_t numRows, size_t numColumns, FILE* stream)
{
  assert(cmr);
  assert(minor);
  assert(stream);

  CMR_CALL( CMRsubmatWriteToStream(cmr, minor->remainingSubmatrix, numRows, numColumns, stream) );

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

  CMR_CALL( CMRminorWriteToStream(cmr, minor, numRows, numColumns, stream) );

  if (stream != stdout)
    fclose(stream);

  return CMR_OKAY;
}


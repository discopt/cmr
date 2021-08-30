#include "regular_internal.h"

#include <cmr/series_parallel.h>

CMR_ERROR CMRregularDecomposeSeriesParallel(CMR* cmr, CMR_DEC** pdec, CMR_CHRMAT** pmatrix, bool ternary)              /**< Whether to consider the signs of the matrix. */
{
  assert(cmr);
  assert(pdec);
  assert(*pdec);
  assert(pmatrix);
  assert(*pmatrix);

//   CMR_CALL( CMRdecomposeBinarySeriesParallel() );

  return CMR_OKAY;
}

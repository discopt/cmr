#include "regular_internal.h"

#include <cmr/graphic.h>
#include <cmr/network.h>

CMR_ERROR CMRregularTestGraphic(CMR* cmr, CMR_CHRMAT** pmatrix, CMR_CHRMAT** ptranspose, bool ternary, bool* pisGraphic,
  CMR_GRAPH** pgraph, CMR_GRAPH_EDGE** pforest, CMR_GRAPH_EDGE** pcoforest, bool** parcsReversed,
  CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(pmatrix);
  assert(ptranspose);
  assert(pisGraphic);

  CMR_CHRMAT* matrix = *pmatrix;
  CMR_CHRMAT* transpose = *ptranspose;

  assert(matrix || transpose);

  if (!transpose)
  {
    CMR_CALL( CMRchrmatTranspose(cmr, matrix, ptranspose) );
    transpose = *ptranspose;
  }

  if (ternary)
  {
    CMR_CALL( CMRtestConetworkMatrix(cmr, transpose, pisGraphic, pgraph, pforest, pcoforest, parcsReversed,
      psubmatrix) );
  }
  else
  {
    CMR_CALL( CMRtestCographicMatrix(cmr, transpose, pisGraphic, pgraph, pforest, pcoforest, psubmatrix) );
  }

  return CMR_OKAY;
}

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "regular_internal.h"
#include "env_internal.h"

bool CMRregularTest(CMR* cmr, CMR_CHRMAT* matrix, int* rowLabels, int* columnLabels,
  CMR_TU_DEC** pdec)
{
  bool isRegular = true;
  bool certify = pdec != NULL;

  assert(cmr);
  assert(matrix);

  /* Perform a 1-sum decomposition. */
  CMR_TU_DEC* dec = NULL;

  int numChildren = CMRregularDecomposeOneSum(cmr, matrix, rowLabels, columnLabels, &dec, true);
  if (certify)
    *pdec = dec;
  if (numChildren <= 1)
  {
    isRegular = CMRregularSequentiallyConnected(cmr, dec, certify, false, false);
  }
  else
  {
    if (certify)
      dec->flags = CMR_TU_DEC_ONE_SUM | CMR_TU_DEC_GRAPHIC | CMR_TU_DEC_COGRAPHIC | CMR_TU_DEC_REGULAR;

    for (int child = 0; child < numChildren; ++child)
    {
      bool result = CMRregularSequentiallyConnected(cmr, dec->children[child], certify,
        false, false);
      isRegular = isRegular && result;
      if (certify)
      {
        dec->flags &= (dec->children[child]->flags
          & (CMR_TU_DEC_REGULAR | CMR_TU_DEC_GRAPHIC | CMR_TU_DEC_COGRAPHIC));
      }
      else if (!result)
        break;
    }
  }

  if (!pdec)
    CMRtudecFree(cmr, &dec);

  return isRegular;
}

bool CMRregularSequentiallyConnected(CMR* cmr, CMR_TU_DEC* decomposition, bool certify, bool notGraphic,
  bool notCographic)
{
  assert(cmr);
  assert(decomposition);
  assert(decomposition->matrix);
  assert(decomposition->transpose);

  return true;
}

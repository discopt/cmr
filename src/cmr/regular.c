#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "regular_internal.h"
#include "env_internal.h"

bool TUregularTest(TU* tu, TU_CHRMAT* matrix, int* rowLabels, int* columnLabels,
  TU_DEC** pdec)
{
  bool isRegular = true;
  bool certify = pdec != NULL;

  assert(tu);
  assert(matrix);

  /* Perform a 1-sum decomposition. */
  TU_DEC* dec = NULL;

  int numChildren = TUregularDecomposeOneSum(tu, matrix, rowLabels, columnLabels, &dec, true);
  if (certify)
    *pdec = dec;
  if (numChildren <= 1)
  {
    isRegular = TUregularSequentiallyConnected(tu, dec, certify, false, false);
  }
  else
  {
    if (certify)
      dec->flags = TU_DEC_ONE_SUM | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC | TU_DEC_REGULAR;

    for (int child = 0; child < numChildren; ++child)
    {
      bool result = TUregularSequentiallyConnected(tu, dec->children[child], certify,
        false, false);
      isRegular = isRegular && result;
      if (certify)
      {
        dec->flags &= (dec->children[child]->flags
          & (TU_DEC_REGULAR | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC));
      }
      else if (!result)
        break;
    }
  }

  if (!pdec)
    TUdecFree(tu, &dec);

  return isRegular;
}

bool TUregularSequentiallyConnected(TU* tu, TU_DEC* decomposition, bool certify, bool notGraphic,
  bool notCographic)
{
  assert(tu);
  assert(decomposition);
  assert(decomposition->matrix);
  assert(decomposition->transpose);

  return true;
}

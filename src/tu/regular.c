#include <tu/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "regular_internal.h"
#include "env_internal.h"

bool TUregularTest(TU* tu, TU_CHRMAT* matrix, int* rowLabels, int* columnLabels,
  TU_DEC** pdecomposition)
{
  bool isRegular = true;
  bool certify = pdecomposition != NULL;

  assert(tu);
  assert(matrix);

  /* Perform a 1-sum decomposition. */
  TU_DEC* decomposition = NULL;

  int numChildren = TUregularDecomposeOneSum(tu, matrix, rowLabels, columnLabels, &decomposition, true);
  if (certify)
    *pdecomposition = decomposition;
  if (numChildren <= 1)
  {
    isRegular = TUregularSequentiallyConnected(tu, decomposition, certify, false, false);
  }
  else
  {
    if (certify)
      decomposition->flags = TU_DEC_ONE_SUM | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC | TU_DEC_REGULAR;

    for (int child = 0; child < numChildren; ++child)
    {
      bool result = TUregularSequentiallyConnected(tu, decomposition->children[child], certify,
        false, false);
      isRegular = isRegular && result;
      if (certify)
      {
        decomposition->flags &= (decomposition->children[child]->flags
          & (TU_DEC_REGULAR | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC));
      }
      else if (!result)
        break;
    }
  }

  if (!pdecomposition)
    TUdecFree(tu, &decomposition);

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

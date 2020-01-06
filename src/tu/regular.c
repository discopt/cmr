#include <tu/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "one_sum.h"

bool TUtestBinaryRegularSequentiallyConnected(TU* tu, TU_CHAR_MATRIX* matrix,
  TU_CHAR_MATRIX* transpose, TU_DEC** decomposition, bool notGraphic, bool notCographic)
{
  return true;
}

int compareComponents(const void* a, const void* b)
{
  return ((TU_ONESUM_COMPONENT*)a)->matrix->numNonzeros -
    ((TU_ONESUM_COMPONENT*)b)->matrix->numNonzeros;
}

bool TUtestBinaryRegular(TU* tu, TU_CHAR_MATRIX* matrix, TU_DEC** decomposition)
{
  bool isRegular = true;
  int numComponents;
  TU_ONESUM_COMPONENT* components;

  assert(tu);
  assert(matrix);

  /* Check entries. */

  assert(TUisBinaryChar(tu, matrix, NULL));

  /* Perform 1-sum decomposition. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  /* Sort components by number of nonzeros. */
  TU_ONESUM_COMPONENT** orderedComponents = 
    (TU_ONESUM_COMPONENT**) malloc( numComponents * sizeof(TU_ONESUM_COMPONENT*) );
  for (int comp = 0; comp < numComponents; ++comp)
    orderedComponents[comp] = &components[comp];
  qsort(orderedComponents, numComponents, sizeof(TU_ONESUM_COMPONENT*), &compareComponents);

  /* Test regularity for each component. */

  if (decomposition)
  {
    *decomposition = (TU_DEC*) malloc( sizeof(TU_DEC) );
  }

  for (int i = 0; i < numComponents; ++i)
  {
    int comp = (orderedComponents[i] - components) / sizeof(TU_ONESUM_COMPONENT*);
    assert(comp >= 0);
    assert(comp < numComponents);

    TU_DEC* compDecomposition = NULL;

    bool compRegular = TUtestBinaryRegularSequentiallyConnected(tu,
      (TU_CHAR_MATRIX*) &components[comp].matrix, (TU_CHAR_MATRIX*) &components[comp].transpose,
      decomposition ? &compDecomposition : NULL, false, false);

    isRegular = isRegular && compRegular;

    if (!decomposition && !compRegular)
      goto cleanup;

    if (decomposition)
    {
      assert(compRegular);
    }
  }

cleanup:
  free(orderedComponents);

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TUfreeCharMatrix(tu, (TU_CHAR_MATRIX**) &components[comp].matrix);
    TUfreeCharMatrix(tu, (TU_CHAR_MATRIX**) &components[comp].transpose);
    TUfreeBlockArray(tu, &components[comp].rowsToOriginal);
    TUfreeBlockArray(tu, &components[comp].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  return isRegular;
}

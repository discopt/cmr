#include <tu/tu.h>

#include "one_sum.h"
#include "sign_internal.h"

#include <stdlib.h>
#include <assert.h>

bool TUtestTotalUnimodularityDouble(TU* tu, TU_SPARSE_DOUBLE* matrix, double epsilon,
  TU_DEC** decomposition, TU_SUBMATRIX** submatrix)
{
  
}

bool TUtestTotalUnimodularityInt(TU* tu, TU_SPARSE_INT* matrix, TU_DEC** decomposition,
  TU_SUBMATRIX** submatrix)
{
  
}

bool TUtestTotalUnimodularityChar(TU* tu, TU_SPARSE_CHAR* matrix, TU_DEC** decomposition,
  TU_SUBMATRIX** submatrix)
{
  int numComponents;
  TU_ONESUM_COMPONENT_CHAR* components;
  bool result;

  assert(tu);
  assert(matrix);

  /* Check entries. */

  if (!TUisTernaryChar(matrix, submatrix))
    return false;

  /* Perform 1-sum decomposition. */

  decomposeOneSumCharToChar(tu, matrix, &numComponents, &components, NULL, NULL, NULL, NULL);

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMATRIX* componentSubmatrix;
    char signFailed = signSequentiallyConnected(tu, &components[comp].matrix,
      &components[comp].transpose, false, submatrix ? &componentSubmatrix : NULL);

    if (signFailed)
    {
      result = false;
      if (submatrix)
        *submatrix = componentSubmatrix;
      goto cleanupComponents;
    }
  }

  /* Check regularity of the binary version. */

cleanupComponents:
  for (int comp = 0; comp < numComponents; ++comp)
  {
    TUclearSparseChar(&components[comp].matrix);
    TUclearSparseChar(&components[comp].transpose);
    free(components[comp].rowsToOriginal);
    free(components[comp].columnsToOriginal);
  }

  return result;
}

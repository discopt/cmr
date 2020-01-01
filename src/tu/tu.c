#include <tu/tu.h>

#include "one_sum.h"
#include "sign_internal.h"

#include <assert.h>

bool TUtestTotalUnimodularityDouble(TU* tu, TU_SPARSE_DOUBLE* matrix, TU_DEC** decomposition,
  TU_SUBMATRIX** submatrix)
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

  assert(tu);
  assert(matrix);

  decomposeOneSumCharToChar(tu, matrix, &numComponents, &components, NULL, NULL, NULL, NULL);

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMATRIX* componentSubmatrix;
    char signFailed = signSequentiallyConnected(tu, &components[comp].matrix,
      &components[comp].transpose, false, submatrix != NULL ? componentSubmatrix : NULL);

    if (signFailed != 0)
    {
      
    }
  }
}

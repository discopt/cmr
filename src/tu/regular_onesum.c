#include <tu/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "decomposition_internal.h"
#include "regular_internal.h"
#include "env_internal.h"
#include "sort.h"
#include "one_sum.h"

int compareOneSumComponents(const void* a, const void* b)
{
  return ((TU_ONESUM_COMPONENT*)a)->matrix->numNonzeros -
    ((TU_ONESUM_COMPONENT*)b)->matrix->numNonzeros;
}

TU_ERROR TUregularDecomposeOneSum(TU* tu, TU_DEC* dec, TU_CHRMAT* matrix)
{
  assert(tu);
  assert(dec);
  assert(matrix);

  /* Perform 1-sum decomposition. */

  size_t numComponents;
  TU_ONESUM_COMPONENT* components = NULL;
  TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  if (numComponents == 1)
  {
    TU_CALL( TUchrmatFree(tu, (TU_CHRMAT**) &components[0].matrix) );
    TU_CALL( TUchrmatFree(tu, (TU_CHRMAT**) &components[0].transpose) );
    TU_CALL( TUfreeBlockArray(tu, &components[0].rowsToOriginal) );
    TU_CALL( TUfreeBlockArray(tu, &components[0].columnsToOriginal) );
  }
  else if (numComponents >= 2)
  {
    /* We create an intermediate array for sorting the components by number of nonzeros. */
    TU_ONESUM_COMPONENT** orderedComponents = NULL;
    TU_CALL( TUallocStackArray(tu, &orderedComponents, numComponents) );
    for (int comp = 0; comp < numComponents; ++comp)
      orderedComponents[comp] = &components[comp];
    TU_CALL( TUsort(tu, numComponents, orderedComponents, sizeof(TU_ONESUM_COMPONENT*), &compareOneSumComponents) );

    /* We now create the children. */
    TU_CALL( TUdecSetNumChildren(tu, dec, numComponents) );
    for (int c = 0; c < numComponents; ++c)
    {
      TU_ONESUM_COMPONENT* component = orderedComponents[c];
      TU_CALL( TUdecCreate(tu, dec, component->matrix->numRows, component->rowsToOriginal,
        component->matrix->numColumns, component->columnsToOriginal, &dec->children[c]) );
      dec->children[c]->matrix = (TU_CHRMAT*) component->matrix;
      dec->children[c]->transpose = (TU_CHRMAT*) component->transpose;
      TU_CALL( TUfreeBlockArray(tu, &component->rowsToOriginal) );
      TU_CALL( TUfreeBlockArray(tu, &component->columnsToOriginal) );
      TU_CALL( TUdecInheritElements(tu, dec->children[c]) );
    }
    dec->type = TU_DEC_ONE_SUM;

    TU_CALL( TUfreeStackArray(tu, &orderedComponents) );
  }

  TU_CALL( TUfreeBlockArray(tu, &components) );

  return TU_OKAY;
}

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "dec_internal.h"
#include "regular_internal.h"
#include "sort.h"
#include "one_sum.h"

int compareOneSumComponents(const void* a, const void* b)
{
  return ((CMR_ONESUM_COMPONENT*)a)->matrix->numNonzeros -
    ((CMR_ONESUM_COMPONENT*)b)->matrix->numNonzeros;
}

CMR_ERROR CMRregularDecomposeOneSum(CMR* cmr, CMR_DEC* dec)
{
  assert(cmr);
  assert(dec);
  assert(dec->matrix);

  /* Perform 1-sum decomposition. */

  size_t numComponents;
  CMR_ONESUM_COMPONENT* components = NULL;
  CMR_CALL( decomposeOneSum(cmr, (CMR_MATRIX*) dec->matrix, sizeof(char), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  if (numComponents == 1)
  {
    CMR_CALL( CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].matrix) );
    CMR_CALL( CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].transpose) );
    CMR_CALL( CMRfreeBlockArray(cmr, &components[0].rowsToOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &components[0].columnsToOriginal) );
  }
  else if (numComponents >= 2)
  {
    /* We create an intermediate array for sorting the components by number of nonzeros. */
    CMR_ONESUM_COMPONENT** orderedComponents = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &orderedComponents, numComponents) );
    for (size_t comp = 0; comp < numComponents; ++comp)
      orderedComponents[comp] = &components[comp];
    CMR_CALL( CMRsort(cmr, numComponents, orderedComponents, sizeof(CMR_ONESUM_COMPONENT*), &compareOneSumComponents) );

    /* We now create the children. */
    CMR_CALL( CMRdecSetNumChildren(cmr, dec, numComponents) );
    for (size_t comp = 0; comp < numComponents; ++comp)
    {
      CMR_ONESUM_COMPONENT* component = orderedComponents[comp];
      CMR_CALL( CMRdecCreate(cmr, dec, component->matrix->numRows, component->rowsToOriginal,
        component->matrix->numColumns, component->columnsToOriginal, &dec->children[comp]) );
      dec->children[comp]->matrix = (CMR_CHRMAT*) component->matrix;
      dec->children[comp]->transpose = (CMR_CHRMAT*) component->transpose;
      CMR_CALL( CMRfreeBlockArray(cmr, &component->rowsToOriginal) );
      CMR_CALL( CMRfreeBlockArray(cmr, &component->columnsToOriginal) );
    }
    dec->type = CMR_DEC_ONE_SUM;

    CMR_CALL( CMRfreeStackArray(cmr, &orderedComponents) );
  }

  CMR_CALL( CMRfreeBlockArray(cmr, &components) );

  return CMR_OKAY;
}

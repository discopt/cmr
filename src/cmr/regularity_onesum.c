#define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "matroid_internal.h"
#include "regularity_internal.h"
#include "sort.h"
#include "one_sum.h"

int compareOneSumComponents(const void* a, const void* b)
{
  return ((CMR_ONESUM_COMPONENT*)b)->matrix->numNonzeros -
    ((CMR_ONESUM_COMPONENT*)a)->matrix->numNonzeros;
}


CMR_ERROR CMRregularitySearchOneSum(CMR* cmr, DecompositionTask* task, DecompositionTask** punprocessed)
{
  assert(cmr);
  assert(task);
  assert(task->dec);
  assert(task->dec->matrix);

#if defined(CMR_DEBUG)
  CMRdbgMsg(2, "Searching for 1-separations.\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, task->dec->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  /* Perform 1-sum decomposition. */

  size_t numComponents;
  CMR_ONESUM_COMPONENT* components = NULL;
  CMR_CALL( decomposeOneSum(cmr, (CMR_MATRIX*) task->dec->matrix, sizeof(char), sizeof(char), &numComponents,
    &components, NULL, NULL, NULL, NULL) );

  if (numComponents == 1)
  {
#if defined(CMR_DEBUG)
    CMRdbgMsg(4, "Matrix is 2-connected.\n", numComponents);
#endif /* CMR_DEBUG */

    CMR_CALL( CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].matrix) );
    CMR_CALL( CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].transpose) );
    CMR_CALL( CMRfreeBlockArray(cmr, &components[0].rowsToOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &components[0].columnsToOriginal) );

    /* Just mark it as 2-connected and add it back to the list of unprocessed tasks. */
    task->dec->testedTwoConnected = true;
    task->next = *punprocessed;
    *punprocessed = task;
  }
  else if (numComponents >= 2)
  {
#if defined(CMR_DEBUG)
    CMRdbgMsg(4, "The 1-sum consists of %zu components.\n", numComponents);
#endif /* CMR_DEBUG */

    /* We create an intermediate array for sorting the components in descending order by number of nonzeros. */

    CMR_ONESUM_COMPONENT** orderedComponents = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &orderedComponents, numComponents) );
    for (size_t comp = 0; comp < numComponents; ++comp)
      orderedComponents[comp] = &components[comp];
    CMR_CALL( CMRsort(cmr, numComponents, orderedComponents, sizeof(CMR_ONESUM_COMPONENT*), &compareOneSumComponents) );

    /* We now create the children. */
    CMR_CALL( CMRmatroiddecUpdateOneSum(cmr, task->dec, numComponents) );
    for (size_t comp = 0; comp < numComponents; ++comp)
    {
      CMR_ONESUM_COMPONENT* component = orderedComponents[comp];
      CMR_CALL( CMRmatroiddecCreateChildFromMatrices(cmr, task->dec, comp, (CMR_CHRMAT*) component->matrix,
        (CMR_CHRMAT*) component->transpose, component->rowsToOriginal, component->columnsToOriginal) );

      CMR_CALL( CMRfreeBlockArray(cmr, &component->rowsToOriginal) );
      CMR_CALL( CMRfreeBlockArray(cmr, &component->columnsToOriginal) );
    }
    task->dec->type = CMR_MATROID_DEC_TYPE_ONE_SUM;

    CMR_CALL( CMRfreeStackArray(cmr, &orderedComponents) );

    for (size_t child = 0; child < numComponents; ++child)
    {
      task->dec->children[child]->testedTwoConnected = true;
      DecompositionTask* childTask = NULL;
      CMR_CALL( CMRregularityTaskCreateRoot(cmr, task->dec->children[child], &childTask, task->params, task->stats,
        task->startClock, task->timeLimit) );
      childTask->next = *punprocessed;
      *punprocessed = childTask;
    }
  }

  CMR_CALL( CMRfreeBlockArray(cmr, &components) );

  return CMR_OKAY;
}

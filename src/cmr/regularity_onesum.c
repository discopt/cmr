// #define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "matroid_internal.h"
#include "regularity_internal.h"
#include "sort.h"
#include "block_decomposition.h"

int compareOneSumComponents(const void* a, const void* b)
{
  return ((CMR_BLOCK*)b)->matrix->numNonzeros -
    ((CMR_BLOCK*)a)->matrix->numNonzeros;
}


CMR_ERROR CMRregularitySearchOneSum(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(task->dec);
  assert(task->dec->matrix);

#if defined(CMR_DEBUG)
  CMRdbgMsg(6, "Searching for 1-separations for the following matrix:\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, task->dec->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  /* Perform 1-sum decomposition. */

  size_t numComponents;
  CMR_BLOCK* components = NULL;
  CMR_CALL( CMRdecomposeBlocks(cmr, (CMR_MATRIX*) task->dec->matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL) );

  if (numComponents == 1)
  {
    CMRdbgMsg(6, "Matrix is 2-connected.\n", numComponents);

    CMR_CALL( CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].matrix) );
    CMR_CALL( CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].transpose) );
    CMR_CALL( CMRfreeBlockArray(cmr, &components[0].rowsToOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &components[0].columnsToOriginal) );

    /* Just mark it as 2-connected and add it back to the list of unprocessed tasks. */
    task->dec->testedTwoConnected = true;
    CMRregularityQueueAdd(queue, task);
  }
  else if (numComponents >= 2)
  {
    CMRdbgMsg(6, "The 1-sum consists of %zu components.\n", numComponents);

    /* We now create the children. */
    CMR_CALL( CMRmatroiddecUpdateOneSum(cmr, task->dec, numComponents) );
    for (size_t comp = 0; comp < numComponents; ++comp)
    {
      CMR_BLOCK* component = &components[comp];
      CMR_CALL( CMRmatroiddecCreateChildFromMatrices(cmr, task->dec, comp, (CMR_CHRMAT*) component->matrix,
        (CMR_CHRMAT*) component->transpose, component->rowsToOriginal, component->columnsToOriginal) );

      CMR_CALL( CMRfreeBlockArray(cmr, &component->rowsToOriginal) );
      CMR_CALL( CMRfreeBlockArray(cmr, &component->columnsToOriginal) );
    }
    task->dec->type = CMR_MATROID_DEC_TYPE_ONE_SUM;

    for (size_t c = numComponents; c > 0; --c)
    {
      size_t child = c - 1;
      task->dec->children[child]->testedTwoConnected = true;
      DecompositionTask* childTask = NULL;
      CMR_CALL( CMRregularityTaskCreateRoot(cmr, task->dec->children[child], &childTask, task->params, task->stats,
        task->startClock, task->timeLimit) );
      CMRregularityQueueAdd(queue, childTask);
    }

    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }

  CMR_CALL( CMRfreeBlockArray(cmr, &components) );

  return CMR_OKAY;
}

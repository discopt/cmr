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
  assert(task->node);
  assert(task->node->matrix);

#if defined(CMR_DEBUG)
  CMRdbgMsg(6, "Searching for 1-separations for the following matrix:\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, task->dec->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  /* Perform 1-sum decomposition. */

  size_t numComponents;
  CMR_BLOCK* components = NULL;
  CMR_CALL( CMRdecomposeBlocks(cmr, (CMR_MATRIX*) task->node->matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL) );

  if (numComponents == 1)
  {
    CMRdbgMsg(6, "Matrix is 2-connected.\n", numComponents);

    CMR_CALL( CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].matrix) );
    CMR_CALL( CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].transpose) );
    CMR_CALL( CMRfreeBlockArray(cmr, &components[0].rowsToOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &components[0].columnsToOriginal) );

    /* Just mark it as 2-connected and add it back to the list of unprocessed tasks. */
    task->node->testedTwoConnected = true;
    CMRregularityQueueAdd(queue, task);
  }
  else if (numComponents >= 2)
  {
    CMRdbgMsg(6, "The 1-sum consists of %zu components.\n", numComponents);

    /* Space for mapping of rows/columns to rows. */
    CMR_ELEMENT* rowsToParent = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &rowsToParent, task->node->numRows) );
    CMR_ELEMENT* columnsToParent = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &columnsToParent, task->node->numColumns) );

    /* We now create the children. */
    CMR_CALL( CMRseymourUpdateOneSum(cmr, task->node, numComponents) );
    for (size_t comp = 0; comp < numComponents; ++comp)
    {
      CMR_BLOCK* component = &components[comp];
      for (size_t row = 0; row < component->matrix->numRows; ++row)
        rowsToParent[row] = CMRrowToElement(component->rowsToOriginal[row]);
      for (size_t column = 0; column < component->matrix->numColumns; ++column)
        columnsToParent[column] = CMRcolumnToElement(component->columnsToOriginal[column]);
      CMR_CALL( CMRseymourCreateChildFromMatrices(cmr, task->node, comp, (CMR_CHRMAT*) component->matrix,
        (CMR_CHRMAT*) component->transpose, rowsToParent, columnsToParent) );

      CMR_CALL( CMRfreeBlockArray(cmr, &component->rowsToOriginal) );
      CMR_CALL( CMRfreeBlockArray(cmr, &component->columnsToOriginal) );
    }
    task->node->type = CMR_SEYMOUR_NODE_TYPE_ONE_SUM;

    CMR_CALL( CMRfreeStackArray(cmr, &columnsToParent) );
    CMR_CALL( CMRfreeStackArray(cmr, &rowsToParent) );

    for (size_t c = numComponents; c > 0; --c)
    {
      size_t child = c - 1;
      task->node->children[child]->testedTwoConnected = true;
      DecompositionTask* childTask = NULL;
      CMR_CALL( CMRregularityTaskCreateRoot(cmr, task->node->children[child], &childTask, task->params, task->stats,
        task->startClock, task->timeLimit) );
      CMRregularityQueueAdd(queue, childTask);
    }

    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }

  CMR_CALL( CMRfreeBlockArray(cmr, &components) );

  return CMR_OKAY;
}

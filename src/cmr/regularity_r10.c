#define CMR_DEBUG /* Uncomment to debug this file. */

#include "env_internal.h"
#include "regularity_internal.h"
#include "matroid_internal.h"

#include <float.h>

CMR_ERROR CMRregularityTestR10(CMR* cmr, DecompositionTask* task, DecompositionTask** punprocessed)
{
  assert(task);
  assert(punprocessed);

  CMR_MATROID_DEC* dec = task->dec;
  assert(dec);

#if defined(CMR_DEBUG)
  CMRdbgMsg(2, "Testing for representation of R10.\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */


  dec->testedR10 = true;
  if (dec->numRows != 5 || dec->numColumns != 5)
    goto cleanup;

  size_t count3 = 0;
  size_t count5 = 0;
  for (size_t row = 0; row < 5; ++row)
  {
    size_t numNonzeros = dec->matrix->rowSlice[row + 1] - dec->matrix->rowSlice[row];  
    if (numNonzeros == 3)
      ++count3;
    else if (numNonzeros == 5)
      ++count5;
    else
      goto cleanup;
  }
  if ((count3 != 4 || count5 != 1) && count3 != 5)
    goto cleanup;

  count3 = 0;
  count5 = 0;
  size_t columnCount[5] = {0, 0, 0, 0, 0};
  for (size_t e = 0; e < dec->matrix->numNonzeros; ++e)
    columnCount[dec->matrix->entryColumns[e]]++;
  for (size_t column = 0; column < 5; ++column)
  {
    if (columnCount[column] == 3)
      ++count3;
    else if (columnCount[column] == 5)
      ++count5;
    else
      goto cleanup;
  }
  if ((count3 != 4 || count5 != 1) && count3 != 5)
    goto cleanup;
  
  /* The number of nonzeros in the rows/columns are 2, 2, 2, 2 and 5. Every 3-connected 5-by-5 matrix with this property
   * represents R10.
   */

  if (dec->isTernary)
  {
    bool isCamion;
    CMR_SUBMAT* violatorSubmatrix = NULL;
    CMR_CALL( CMRtestCamionSigned(cmr, dec->matrix, &isCamion, &violatorSubmatrix,
      task->stats ? &task->stats->network.camion : NULL, DBL_MAX) );

    if (violatorSubmatrix)
    {
      assert(!isCamion);

      /* TODO: Do we actually need the child node or just the irregularity information? */
      CMR_CALL( CMRmatroiddecUpdateSubmatrix(cmr, dec, violatorSubmatrix, CMR_MATROID_DEC_TYPE_DETERMINANT) );
      CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );

      /* Task is done. */
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );

      goto cleanup;
    }
    else
    {
      assert(isCamion);
    }
  }

  dec->type = CMR_MATROID_DEC_TYPE_R10;

cleanup:

  if (dec->type == CMR_MATROID_DEC_TYPE_R10)
  {
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }
  else
  {
    task->next = *punprocessed;
    *punprocessed = task;
  }

  return CMR_OKAY;
}

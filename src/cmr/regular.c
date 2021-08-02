#define CMR_DEBUG /** Uncomment to debug the regularity check. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "decomposition_internal.h"
#include "env_internal.h"
#include "regular_internal.h"

CMR_ERROR CMRtestBinaryRegularConnected(CMR* cmr, CMR_TU_DEC* dec, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose,
  bool checkPlanarity, bool certify, bool constructDecomposition, bool* pisRegular)
{
  assert(cmr);
  assert(dec);
  assert(matrix);

  CMRdbgMsg(2, "Testing binary %dx%d 1-connected matrix for regularity.\n", matrix->numRows, matrix->numColumns);

  

  assert(false);

  return CMR_OKAY;
}

CMR_ERROR CMRtestBinaryRegular(CMR* cmr, CMR_CHRMAT* matrix, CMR_ELEMENT* rowElements, CMR_ELEMENT* columnElements,
  bool checkPlanarity, bool certify, bool *pisRegular, CMR_TU_DEC** pdec)
{
  assert(cmr);
  assert(matrix);

  CMRdbgMsg(0, "Testing binary %dx%d matrix for regularity.\n", matrix->numRows, matrix->numColumns);

  CMR_TU_DEC* dec = NULL;
  CMR_CALL( CMRtudecCreate(cmr, NULL, matrix->numRows, NULL, matrix->numColumns, NULL, &dec) );
  assert(dec);

  /* Either copy or create canonical row elements. */
  if (rowElements)
    CMR_CALL( CMRduplicateBlockArray(cmr, &dec->rowElements, dec->numRows, rowElements) );
  else
  {
    CMR_CALL( CMRallocBlockArray(cmr, &dec->rowElements, dec->numRows) );
    for (size_t row = 0; row < dec->numRows; ++row)
      dec->rowElements[row] = CMRrowToElement(row);
  }

  /* Either copy or create canonical column elements. */
  if (columnElements)
    CMR_CALL( CMRduplicateBlockArray(cmr, &dec->columnElements, dec->numColumns, columnElements) );
  else
  {
    CMR_CALL( CMRallocBlockArray(cmr, &dec->columnElements, dec->numColumns) );
    for (size_t column = 0; column < dec->numColumns; ++column)
      dec->columnElements[column] = CMRcolumnToElement(column);
  }

  CMR_CALL( CMRregularDecomposeOneSum(cmr, dec, matrix) );

  CMRdbgMsg(2, "1-sum decomposition yields %d components.\n", dec->numChildren == 0 ? 1 : dec->numChildren);
#if defined(CMR_DEBUG)
  CMR_CALL( CMRtudecPrint(stdout, dec, 2) );
#endif /* CMR_DEBUG */

  bool isRegular = true;
  if (dec->numChildren)
  {
    for (size_t c = 0; c < dec->numChildren; ++c)
    {
      if (isRegular || pdec)
      {
        bool childIsRegular = true;
        CMR_CALL( CMRtestBinaryRegularConnected(cmr, dec, dec->children[c]->matrix, dec->children[c]->transpose,
          checkPlanarity, certify, pdec, &childIsRegular) );

        if (!childIsRegular)
          isRegular = false;
      }
      else
        CMR_CALL( CMRtudecFree(cmr, &dec->children[c]) );
    }
  }
  else
  {
    CMR_CALL( CMRtestBinaryRegularConnected(cmr, dec, matrix, NULL, checkPlanarity, certify, pdec, &isRegular) );
  }

  if (pisRegular)
    *pisRegular = isRegular;

  CMR_CALL( CMRtudecComputeRegularity(dec) );

  *pdec = dec;

  return CMR_OKAY;
}

// bool TUregularTest(CMR* cmr, CMR_CHRMAT* matrix, int* rowLabels, int* columnLabels,
//   CMR_TU_DEC** pdec)
// {
//   bool isRegular = true;
//   bool certify = pdec != NULL;
// 
//   assert(cmr);
//   assert(matrix);
// 
//   /* Perform a 1-sum decomposition. */
//   CMR_TU_DEC* dec = NULL;
// 
//   int numChildren = TUregularDecomposeOneSum(cmr, matrix, rowLabels, columnLabels, &dec, true);
//   if (certify)
//     *pdec = dec;
//   if (numChildren <= 1)
//   {
//     isRegular = TUregularSequentiallyConnected(cmr, dec, certify, false, false);
//   }
//   else
//   {
//     if (certify)
//       dec->flags = CMR_TU_DEC_ONE_SUM | CMR_TU_DEC_GRAPHIC | CMR_TU_DEC_COGRAPHIC | CMR_TU_DEC_REGULAR;
// 
//     for (int child = 0; child < numChildren; ++child)
//     {
//       bool result = TUregularSequentiallyConnected(cmr, dec->children[child], certify,
//         false, false);
//       isRegular = isRegular && result;
//       if (certify)
//       {
//         dec->flags &= (dec->children[child]->flags
//           & (CMR_TU_DEC_REGULAR | CMR_TU_DEC_GRAPHIC | CMR_TU_DEC_COGRAPHIC));
//       }
//       else if (!result)
//         break;
//     }
//   }
// 
//   if (!pdec)
//     CMRtudecFree(cmr, &dec);
// 
//   return isRegular;
// }
// 
// bool TUregularSequentiallyConnected(CMR* cmr, CMR_TU_DEC* decomposition, bool certify, bool notGraphic,
//   bool notCographic)
// {
//   assert(cmr);
//   assert(decomposition);
//   assert(decomposition->matrix);
//   assert(decomposition->transpose);
// 
//   return true;
// }

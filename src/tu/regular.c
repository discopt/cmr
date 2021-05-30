#define TU_DEBUG /** Uncomment to debug the regularity check. */

#include <tu/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "decomposition_internal.h"
#include "env_internal.h"
#include "regular_internal.h"

TU_ERROR TUtestBinaryRegularConnected(TU* tu, TU_DEC* dec, TU_CHRMAT* matrix, TU_CHRMAT* transpose,
  bool checkPlanarity, bool certify, bool constructDecomposition, bool* pisRegular)
{
  assert(tu);
  assert(dec);
  assert(matrix);

  TUdbgMsg(2, "Testing binary %dx%d 1-connected matrix for regularity.\n", matrix->numRows, matrix->numColumns);

  

  assert(false);

  return TU_OKAY;
}

TU_ERROR TUtestBinaryRegular(TU* tu, TU_CHRMAT* matrix, TU_ELEMENT* rowElements, TU_ELEMENT* columnElements,
  bool checkPlanarity, bool certify, bool *pisRegular, TU_DEC** pdec)
{
  assert(tu);
  assert(matrix);

  TUdbgMsg(0, "Testing binary %dx%d matrix for regularity.\n", matrix->numRows, matrix->numColumns);

  TU_DEC* dec = NULL;
  TU_CALL( TUdecCreate(tu, NULL, matrix->numRows, NULL, matrix->numColumns, NULL, &dec) );
  assert(dec);

  /* Either copy or create canonical row elements. */
  if (rowElements)
    TU_CALL( TUduplicateBlockArray(tu, &dec->rowElements, dec->numRows, rowElements) );
  else
  {
    TU_CALL( TUallocBlockArray(tu, &dec->rowElements, dec->numRows) );
    for (size_t row = 0; row < dec->numRows; ++row)
      dec->rowElements[row] = TUrowToElement(row);
  }

  /* Either copy or create canonical column elements. */
  if (columnElements)
    TU_CALL( TUduplicateBlockArray(tu, &dec->columnElements, dec->numColumns, columnElements) );
  else
  {
    TU_CALL( TUallocBlockArray(tu, &dec->columnElements, dec->numColumns) );
    for (size_t column = 0; column < dec->numColumns; ++column)
      dec->columnElements[column] = TUcolumnToElement(column);
  }

  TU_CALL( TUregularDecomposeOneSum(tu, dec, matrix) );

  TUdbgMsg(2, "1-sum decomposition yields %d components.\n", dec->numChildren == 0 ? 1 : dec->numChildren);
#if defined(TU_DEBUG)
  TU_CALL( TUdecPrint(stdout, dec, 2) );
#endif /* TU_DEBUG */

  bool isRegular = true;
  if (dec->numChildren)
  {
    for (size_t c = 0; c < dec->numChildren; ++c)
    {
      if (isRegular || pdec)
      {
        bool childIsRegular = true;
        TU_CALL( TUtestBinaryRegularConnected(tu, dec, dec->children[c]->matrix, dec->children[c]->transpose,
          checkPlanarity, certify, pdec, &childIsRegular) );

        if (!childIsRegular)
          isRegular = false;
      }
      else
        TU_CALL( TUdecFree(tu, &dec->children[c]) );
    }
  }
  else
  {
    TU_CALL( TUtestBinaryRegularConnected(tu, dec, matrix, NULL, checkPlanarity, certify, pdec, &isRegular) );
  }

  if (pisRegular)
    *pisRegular = isRegular;

  TU_CALL( TUdecComputeRegularity(dec) );

  *pdec = dec;

  return TU_OKAY;
}

// bool TUregularTest(TU* tu, TU_CHRMAT* matrix, int* rowLabels, int* columnLabels,
//   TU_DEC** pdec)
// {
//   bool isRegular = true;
//   bool certify = pdec != NULL;
// 
//   assert(tu);
//   assert(matrix);
// 
//   /* Perform a 1-sum decomposition. */
//   TU_DEC* dec = NULL;
// 
//   int numChildren = TUregularDecomposeOneSum(tu, matrix, rowLabels, columnLabels, &dec, true);
//   if (certify)
//     *pdec = dec;
//   if (numChildren <= 1)
//   {
//     isRegular = TUregularSequentiallyConnected(tu, dec, certify, false, false);
//   }
//   else
//   {
//     if (certify)
//       dec->flags = TU_DEC_ONE_SUM | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC | TU_DEC_REGULAR;
// 
//     for (int child = 0; child < numChildren; ++child)
//     {
//       bool result = TUregularSequentiallyConnected(tu, dec->children[child], certify,
//         false, false);
//       isRegular = isRegular && result;
//       if (certify)
//       {
//         dec->flags &= (dec->children[child]->flags
//           & (TU_DEC_REGULAR | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC));
//       }
//       else if (!result)
//         break;
//     }
//   }
// 
//   if (!pdec)
//     TUdecFree(tu, &dec);
// 
//   return isRegular;
// }
// 
// bool TUregularSequentiallyConnected(TU* tu, TU_DEC* decomposition, bool certify, bool notGraphic,
//   bool notCographic)
// {
//   assert(tu);
//   assert(decomposition);
//   assert(decomposition->matrix);
//   assert(decomposition->transpose);
// 
//   return true;
// }

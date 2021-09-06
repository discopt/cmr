#define CMR_DEBUG /* Uncomment to debug this file. */

#include "dec_internal.h"
#include "regular_internal.h"
#include "env_internal.h"

// CMR_ERROR CMRregularListMatrixFree(CMR* cmr, ListMatrix** pmatrix)
// {
//   assert(cmr);
//   assert(pmatrix);
// 
//   ListMatrix* matrix = *pmatrix;
//   if (!matrix)
//     return CMR_OKAY;
// 
//   CMR_CALL( CMRfreeBlockArray(cmr, &matrix->columns) );
//   CMR_CALL( CMRfreeBlockArray(cmr, &matrix->rows) );
//   CMR_CALL( CMRfreeBlock(cmr, pmatrix) );
// 
//   return CMR_OKAY;
// }

CMR_ERROR CMRregularConstructNestedMinorSequence(CMR* cmr, CMR_DEC* dec, bool ternary, CMR_SUBMAT* wheelSubmatrix,
//   ListMatrix** pnestedMatrix,
  NestedMinor** pnestedMinors, size_t* pnumNestedMinors, CMR_SUBMAT** psubmatrix,
  CMR_REGULAR_PARAMETERS* params)
{
  assert(cmr);
  assert(dec);
  assert(wheelSubmatrix);
//   assert(pnestedMatrix);
//   assert(pnestedMinors);
  assert(pnumNestedMinors);
  assert(params);

  CMRdbgMsg(4, "Attempting to construct a sequence of 3-connected nested minors.\n");


  /* Create list representation of the wheel submatrix. */
//   CMR_CALL( CMRallocBlock(cmr, pnestedMatrix) );
//   ListMatrix* nestedMatrix = *pnestedMatrix;

  /* Mapping of rows to original rows. */
//   nestedMatrix->numRows = wheelSubmatrix->numRows;
//   nestedMatrix->rows = NULL;
//   CMR_CALL( CMRallocBlockArray(cmr, &nestedMatrix->rows, dec->matrix->numRows) );
//   for (size_t row = 0; row < wheelSubmatrix->numRows; ++row)
//   {
//     nestedMatrix->rows[row].element = CMRrowToElement(wheelSubmatrix->rows[row]);
//     nestedMatrix->rows[row].numNonzeros = 0;
//   }

  /* Mapping of columns to original columns. */
//   nestedMatrix->numColumns = wheelSubmatrix->numColumns;
//   nestedMatrix->columns = NULL;
//   CMR_CALL( CMRallocBlockArray(cmr, &nestedMatrix->columns, dec->matrix->numColumns) );
//   for (size_t column = 0; column < wheelSubmatrix->numColumns; ++column)
//   {
//     nestedMatrix->columns[column].element = CMRcolumnToElement(wheelSubmatrix->columns[column]);
//     nestedMatrix->columns[column].numNonzeros = 0;
//   }


  assert("Construction of nested 3-connected minor sequence not fully implemented." == 0);

  return CMR_OKAY;
}

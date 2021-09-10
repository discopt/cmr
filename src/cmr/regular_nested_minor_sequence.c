#define CMR_DEBUG /* Uncomment to debug this file. */

#include "dec_internal.h"
#include "regular_internal.h"
#include "env_internal.h"

#include "hashtable.h"
#include "listmatrix.h"

typedef struct
{
  long long hashValue;                /**< \brief Hash value of this element. */
  CMR_LISTHASHTABLE_ENTRY hashEntry;  /**< \brief Entry in row or column hashtable. */
} ElementData;


/**
 * \brief Returns the smallest power of 2 at least as large as \p x.
 */

static
size_t nextPower2(size_t x)
{
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return x + 1;
}

static
CMR_ERROR createHashVector(
  CMR* cmr,                 /**< \ref CMR environment. */
  long long** phashVector,  /**< Pointer for storing the hash vector. */
  size_t size               /**< Size of hash vector. */
)
{
  assert(cmr);

  CMR_CALL( CMRallocStackArray(cmr, phashVector, size) );
  long long* hashVector = *phashVector;
  size_t h = 1;
  for (size_t e = 0; e < size; ++e)
  {
    hashVector[e] = h;
    CMRdbgMsg(2, "Entry %d has hash %ld.\n", e, h);
    h = projectSignedHash(3 * h);
  }

  return CMR_OKAY;
}

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

  size_t numRows = dec->matrix->numRows;
  size_t numColumns = dec->matrix->numColumns;

  CMRdbgMsg(4, "Attempting to construct a sequence of 3-connected nested minors.\n");

  ListMatrix* nestedMatrix = NULL;
  CMR_CALL( CMRlistmatrixAlloc(cmr, numRows, numColumns, 2 * dec->matrix->numNonzeros, &nestedMatrix) );
  CMR_CALL( CMRlistmatrixInitializeFromSubmatrix(cmr, nestedMatrix, dec->matrix, wheelSubmatrix) );

  ListMatrix* remainingMatrix = NULL;
  CMR_CALL( CMRlistmatrixAlloc(cmr, numRows, numColumns, 2 * dec->matrix->numNonzeros, &remainingMatrix) );
  CMR_CALL( CMRlistmatrixInitializeFromSubmatrixComplement(cmr, remainingMatrix, dec->matrix, wheelSubmatrix) );

  CMR_LISTHASHTABLE* rowHashtable = NULL;
  CMR_CALL( CMRlisthashtableCreate(cmr, &rowHashtable, nextPower2(numRows), numRows ) );

  CMR_LISTHASHTABLE* columnHashtable = NULL;
  CMR_CALL( CMRlisthashtableCreate(cmr, &columnHashtable, nextPower2(numColumns), numColumns ) );

  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector, numRows > numColumns ? numRows : numColumns) );

  ElementData* rowData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowData, numRows) );
  for (size_t row = 0; row < numRows; ++row)
  {
    rowData[row].hashValue = 0;
    rowData[row].hashEntry = SIZE_MAX;
  }

  ElementData* columnData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnData, dec->matrix->numColumns) );
  for (size_t column = 0; column < numColumns; ++column)
  {
    columnData[column].hashValue = 0;
    columnData[column].hashEntry = SIZE_MAX;
  }

  do
  {
    printf("Nested matrix:\n");
    CMR_CALL( CMRlistmatrixPrintDense(cmr, nestedMatrix, stdout) );

    printf("Remaining matrix:\n");
    CMR_CALL( CMRlistmatrixPrintDense(cmr, remainingMatrix, stdout) );


  }
  while (false);

  CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowData) );
  CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );

  CMR_CALL( CMRlisthashtableFree(cmr, &columnHashtable) );
  CMR_CALL( CMRlisthashtableFree(cmr, &rowHashtable) );

  CMR_CALL( CMRlistmatrixFree(cmr, &remainingMatrix) );
  CMR_CALL( CMRlistmatrixFree(cmr, &nestedMatrix) );

  assert("Construction of nested 3-connected minor sequence not fully implemented." == 0);

  return CMR_OKAY;
}

#define CMR_DEBUG /* Uncomment to debug this file. */

#include "env_internal.h"
#include "regular_internal.h"
#include "dec_internal.h"

CMR_ERROR CMRregularThreeConnectedIsR10(CMR* cmr, CMR_DEC* dec, bool* pisR10)
{
  assert(cmr);
  assert(dec);
  assert(pisR10);

  *pisR10 = false;
  if (dec->matrix->numRows != 5 || dec->matrix->numColumns != 5)
    return CMR_OKAY;

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
      return CMR_OKAY;
  }
  if ((count3 != 4 || count5 != 1) && count3 != 5)
    return CMR_OKAY;

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
      return CMR_OKAY;
  }
  if ((count3 != 4 || count5 != 1) && count3 != 5)
    return CMR_OKAY;
  
  /* The number of nonzeros in the rows/columns are 2, 2, 2, 2 and 5. Every 3-connected 5-by-5 matrix with this property
   * represents R10.
   */

  dec->type = CMR_DEC_SPECIAL_R10;
  *pisR10 = true;

  return CMR_OKAY;
}

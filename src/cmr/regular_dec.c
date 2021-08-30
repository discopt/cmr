#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"

void CMRcreateDec(CMR* cmr, CMR_DEC** pdec)
{
  assert(pdec != NULL);
  assert(*pdec == NULL);

  CMRallocBlock(cmr, pdec);
  CMR_DEC* dec = *pdec;
  dec->matrix = NULL;
  dec->transpose = NULL;
  dec->rowLabels = NULL;
  dec->columnLabels = NULL;
  dec->parentRows = NULL;
  dec->parentColumns = NULL;
  dec->flags = 0;
  dec->children = NULL;
  dec->numChildren = 0;
  dec->graph = NULL;
  dec->cograph = NULL;
}

void CMRdecFree(CMR* cmr, CMR_DEC** pdec)
{
  assert(pdec);
  assert(*pdec);

  CMR_DEC* dec = *pdec;

  if (dec->numChildren > 0)
  {
    for (int c = 0; c < dec->numChildren; ++c)
      CMRdecFree(cmr, &dec->children[c]);
    CMRfreeBlockArray(cmr, &dec->children);
  }
  if (dec->matrix)
    CMRchrmatFree(cmr, &dec->matrix);
  if (dec->transpose)
    CMRchrmatFree(cmr, &dec->transpose);
  if (dec->rowLabels)
    CMRfreeBlockArray(cmr, &dec->rowLabels);
  if (dec->columnLabels)
    CMRfreeBlockArray(cmr, &dec->columnLabels);
  if (dec->parentRows)
    CMRfreeBlockArray(cmr, &dec->parentRows);
  if (dec->parentColumns)
    CMRfreeBlockArray(cmr, &dec->parentColumns);
  if (dec->graph)
    CMRgraphFree(cmr, &dec->graph);
  if (dec->cograph)
    CMRgraphFree(cmr, &dec->cograph);

  CMRfreeBlock(cmr, pdec);
}

bool CMRdecIsLeaf(CMR_DEC* dec)
{
  assert(dec);
  return dec->numChildren == 0;
}

bool CMRdecIsRegular(CMR_DEC* dec)
{
  assert(dec);
  return dec->flags & CMR_DEC_REGULAR;
}

bool CMRdecIsGraphic(CMR_DEC* dec)
{
  assert(dec);
  return dec->flags & CMR_DEC_GRAPHIC;
}

bool CMRdecIsCographic(CMR_DEC* dec)
{
  assert(dec);
  return dec->flags & CMR_DEC_COGRAPHIC;
}

char CMRdecIsSum(CMR_DEC* dec)
{
  assert(dec);
  char result = dec->flags & CMR_DEC_TYPE_MASK;
  return result <= 3 ? result : 0;
}

int CMRdecNumRows(CMR_DEC* dec)
{
  assert(dec);
  assert(dec->matrix);
  return dec->matrix->numRows;
}

int CMRdecNumColumns(CMR_DEC* dec)
{
  assert(dec);
  assert(dec->matrix);
  return dec->matrix->numColumns;
}

int CMRgetDecRankLowerLeft(CMR_DEC* dec)
{
  assert(dec);

  if ((dec->flags & CMR_DEC_TYPE_MASK) == CMR_DEC_THREE_SUM
    && !(dec->flags & CMR_DEC_RANK_UPPER_RIGHT))
    return 2;
  else if (dec->flags & (CMR_DEC_RANK_LOWER_LEFT))
    return 1;
  else
    return 0;
}

int CMRgetDecRankUpperRight(CMR_DEC* dec)
{
  assert(dec);

  if ((dec->flags & CMR_DEC_TYPE_MASK) == CMR_DEC_THREE_SUM
    && !(dec->flags & CMR_DEC_RANK_LOWER_LEFT))
    return 2;
  else if (dec->flags & (CMR_DEC_RANK_UPPER_RIGHT))
    return 1;
  else
    return 0;
}

void CMRcreateDecChild(CMR* cmr, CMR_DEC* dec, int numRows, int* rows, int numColumns, int* columns,
  int numNonzeros, int numExtraRows, int numExtraColumns, int numExtraNonzeros,
  bool constructDecomposition, CMR_DEC** presult)
{
  assert(cmr);
  assert(dec);
  assert(rows);
  assert(columns);

  CMR_CHRMAT* parentMatrix = dec->matrix;
  assert(parentMatrix);

  CMRcreateDec(cmr, presult);
  CMR_DEC* result = *presult;

  /* Copy parentRows from rows. */

  CMRallocBlockArray(cmr, &result->parentRows, numRows + numExtraRows);
  for (int row = 0; row < numRows; ++row)
    result->parentRows[row] = rows[row];
  for (int row = 0; row < numExtraRows; ++row)
    result->parentRows[numRows + row] = 0;

  /* Copy parentColumns from columns. */

  CMRallocBlockArray(cmr, &result->parentColumns, numColumns + numExtraColumns);
  for (int column = 0; column < numColumns; ++column)
    result->parentColumns[column] = columns[column];
  for (int column = 0; column < numExtraColumns; ++column)
    result->parentColumns[numColumns + column] = 0;

  /* Create the child matrix. */

  CMRchrmatCreate(cmr, &result->matrix, numRows, numColumns, 0);
  CMR_CHRMAT* childMatrix = result->matrix;

  int* columnMap = NULL;
  CMRallocStackArray(cmr, &columnMap, parentMatrix->numColumns);
  for (int column = 0; column < parentMatrix->numColumns; ++column)
    columnMap[column] = -1;
  for (int i = 0; i < numColumns; ++i)
    columnMap[columns[i]] = i;

  /* If not provided, we count nonzeros in submatrix of parent matrix. */

  if (numNonzeros == 0)
  {
    for (int i = 0; i < numRows; ++i)
    {
      size_t parentRow = rows[i];
      size_t first = parentMatrix->rowSlice[parentRow];
      size_t beyond = parentMatrix->rowSlice[parentRow + 1];
      for (size_t entry = first; entry < beyond; ++entry)
      {
        size_t parentColumn = parentMatrix->entryColumns[entry];
        if (columnMap[parentColumn] >= 0)
          ++numNonzeros;
      }
    }
  }

  CMRallocBlockArray(cmr, &childMatrix->entryColumns, numNonzeros + numExtraNonzeros);
  CMRallocBlockArray(cmr, &childMatrix->entryValues, numNonzeros + numExtraNonzeros);

  /* Write nonzeros to child matrix. */

  numNonzeros = 0;
  size_t countNonzeros = 0; // TODO: Added to make it compile.
  for (size_t childRow = 0; childRow < numRows; ++childRow)
  {
    size_t parentRow = rows[childRow];
    size_t first = parentMatrix->rowSlice[parentRow];
    size_t beyond = parentMatrix->rowSlice[parentRow + 1];
    childMatrix->rowSlice[childRow] = countNonzeros;
    for (size_t entry = first; entry < beyond; ++entry)
    {
      size_t parentColumn = parentMatrix->entryColumns[entry];
      size_t childColumn = columnMap[parentColumn];
      if (childColumn >= 0)
      {
        childMatrix->entryColumns[countNonzeros] = childColumn;
        childMatrix->entryValues[countNonzeros] = parentMatrix->entryValues[entry];
        ++childColumn;
      }
    }
  }

  CMRfreeStackArray(cmr, &columnMap);

  if (constructDecomposition)
  {
    CMRchrmatTranspose(cmr, result->matrix, &result->transpose);
  }
}

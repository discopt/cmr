#include <tu/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"

void TUcreateDec(TU* tu, TU_DEC** pdec)
{
  assert(pdec != NULL);
  assert(*pdec == NULL);

  TUallocBlock(tu, pdec);
  TU_DEC* dec = *pdec;
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

void TUdecFree(TU* tu, TU_DEC** pdec)
{
  assert(pdec);
  assert(*pdec);

  TU_DEC* dec = *pdec;

  if (dec->numChildren > 0)
  {
    for (int c = 0; c < dec->numChildren; ++c)
      TUdecFree(tu, &dec->children[c]);
    TUfreeBlockArray(tu, &dec->children);
  }
  if (dec->matrix)
    TUchrmatFree(tu, &dec->matrix);
  if (dec->transpose)
    TUchrmatFree(tu, &dec->transpose);
  if (dec->rowLabels)
    TUfreeBlockArray(tu, &dec->rowLabels);
  if (dec->columnLabels)
    TUfreeBlockArray(tu, &dec->columnLabels);
  if (dec->parentRows)
    TUfreeBlockArray(tu, &dec->parentRows);
  if (dec->parentColumns)
    TUfreeBlockArray(tu, &dec->parentColumns);
  if (dec->graph)
    TUgraphFree(tu, &dec->graph);
  if (dec->cograph)
    TUgraphFree(tu, &dec->cograph);

  TUfreeBlock(tu, pdec);
}

bool TUdecIsLeaf(TU_DEC* dec)
{
  assert(dec);
  return dec->numChildren == 0;
}

bool TUdecIsRegular(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_REGULAR;
}

bool TUdecIsGraphic(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_GRAPHIC;
}

bool TUdecIsCographic(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_COGRAPHIC;
}

char TUdecIsSum(TU_DEC* dec)
{
  assert(dec);
  char result = dec->flags & TU_DEC_TYPE_MASK;
  return result <= 3 ? result : 0;
}

int TUdecNumRows(TU_DEC* dec)
{
  assert(dec);
  assert(dec->matrix);
  return dec->matrix->numRows;
}

int TUdecNumColumns(TU_DEC* dec)
{
  assert(dec);
  assert(dec->matrix);
  return dec->matrix->numColumns;
}

int TUgetDecRankLowerLeft(TU_DEC* dec)
{
  assert(dec);

  if ((dec->flags & TU_DEC_TYPE_MASK) == TU_DEC_THREE_SUM 
    && !(dec->flags & TU_DEC_RANK_UPPER_RIGHT))
    return 2;
  else if (dec->flags & (TU_DEC_RANK_LOWER_LEFT))
    return 1;
  else
    return 0;
}

int TUgetDecRankUpperRight(TU_DEC* dec)
{
  assert(dec);

  if ((dec->flags & TU_DEC_TYPE_MASK) == TU_DEC_THREE_SUM 
    && !(dec->flags & TU_DEC_RANK_LOWER_LEFT))
    return 2;
  else if (dec->flags & (TU_DEC_RANK_UPPER_RIGHT))
    return 1;
  else
    return 0;
}

void TUcreateDecChild(TU* tu, TU_DEC* dec, int numRows, int* rows, int numColumns, int* columns,
  int numNonzeros, int numExtraRows, int numExtraColumns, int numExtraNonzeros,
  bool constructDecomposition, TU_DEC** presult)
{
  assert(tu);
  assert(dec);
  assert(rows);
  assert(columns);

  TU_CHRMAT* parentMatrix = dec->matrix;
  assert(parentMatrix);

  TUcreateDec(tu, presult);
  TU_DEC* result = *presult;

  /* Copy parentRows from rows. */

  TUallocBlockArray(tu, &result->parentRows, numRows + numExtraRows);
  for (int row = 0; row < numRows; ++row)
    result->parentRows[row] = rows[row];
  for (int row = 0; row < numExtraRows; ++row)
    result->parentRows[numRows + row] = 0;

  /* Copy parentColumns from columns. */

  TUallocBlockArray(tu, &result->parentColumns, numColumns + numExtraColumns);
  for (int column = 0; column < numColumns; ++column)
    result->parentColumns[column] = columns[column];
  for (int column = 0; column < numExtraColumns; ++column)
    result->parentColumns[numColumns + column] = 0;

  /* Create the child matrix. */
  
  TUchrmatCreate(tu, &result->matrix, numRows, numColumns, 0);
  TU_CHRMAT* childMatrix = result->matrix;

  int* columnMap = NULL;
  TUallocStackArray(tu, &columnMap, parentMatrix->numColumns);
  for (int column = 0; column < parentMatrix->numColumns; ++column)
    columnMap[column] = -1;
  for (int i = 0; i < numColumns; ++i)
    columnMap[columns[i]] = i;

  /* If not provided, we count nonzeros in submatrix of parent matrix. */

  if (numNonzeros == 0)
  {
    for (int i = 0; i < numRows; ++i)
    {
      int parentRow = rows[i];
      int begin = parentMatrix->rowStarts[parentRow];
      int end = parentRow + 1 < parentMatrix->numRows ? parentMatrix->rowStarts[parentRow + 1]
        : parentMatrix->numNonzeros;
      for (int entry = begin; entry < end; ++entry)
      {
        int parentColumn = parentMatrix->entryColumns[entry];
        if (columnMap[parentColumn] >= 0)
          ++numNonzeros;
      }
    }
  }

  TUallocBlockArray(tu, &childMatrix->entryColumns, numNonzeros + numExtraNonzeros);
  TUallocBlockArray(tu, &childMatrix->entryValues, numNonzeros + numExtraNonzeros);

  /* Write nonzeros to child matrix. */

  numNonzeros = 0;
  int countNonzeros = 0; // TODO: Added to make it compile.
  for (int childRow = 0; childRow < numRows; ++childRow)
  {
    int parentRow = rows[childRow];
    int begin = parentMatrix->rowStarts[parentRow];
    int end = parentRow + 1 < parentMatrix->numRows ? parentMatrix->rowStarts[parentRow + 1]
      : parentMatrix->numNonzeros;
    childMatrix->rowStarts[childRow] = countNonzeros;
    for (int entry = begin; entry < end; ++entry)
    {
      int parentColumn = parentMatrix->entryColumns[entry];
      int childColumn = columnMap[parentColumn];
      if (childColumn >= 0)
      {
        childMatrix->entryColumns[countNonzeros] = childColumn;
        childMatrix->entryValues[countNonzeros] = parentMatrix->entryValues[entry];
        ++childColumn;
      }
    }
  }

  TUfreeStackArray(tu, &columnMap);

  if (constructDecomposition)
  {
    TUchrmatTranspose(tu, result->matrix, &result->transpose);
  }
}

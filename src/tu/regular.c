#include <tu/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "one_sum.h"


void TUfreeDec(TU* tu, TU_DEC** dec)
{
  assert(dec);
  assert(*dec);
  if ((*dec)->numChildren > 0)
  {
    for (int c = 0; c < (*dec)->numChildren; ++c)
      TUfreeDec(tu, &(*dec)->children[c]);
    TUfreeBlockArray(tu, &(*dec)->children);
  }
  if ((*dec)->matrix)
    TUfreeCharMatrix(tu, &(*dec)->matrix);
  if ((*dec)->transpose)
    TUfreeCharMatrix(tu, &(*dec)->transpose);
  if ((*dec)->graph)
    TUfreeGraph(tu, &(*dec)->graph);
  if ((*dec)->cograph)
    TUfreeGraph(tu, &(*dec)->cograph);
  TUfreeBlock(tu, dec);
}

bool TUisDecLeaf(TU_DEC* dec)
{
  assert(dec);
  return dec->numChildren == 0;
}

bool TUisDecRegular(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_REGULAR;
}

bool TUisDecGraphic(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_GRAPHIC;
}

bool TUisDecCographic(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_COGRAPHIC;
}

char TUisDecSum(TU_DEC* dec)
{
  assert(dec);
  char result = dec->flags & TU_DEC_TYPE_MASK;
  return result <= 3 ? result : 0;
}

int TUgetDecNumRows(TU_DEC* dec)
{
  assert(dec);
  assert(dec->matrix);
  return dec->matrix->numRows;
}

int TUgetDecNumColumns(TU_DEC* dec)
{
  assert(dec);
  assert(dec->matrix);
  return dec->matrix->numColumns;
}

/**
 * \brief Returns the rank of the lower-left submatrix of this node of the tree.
 *
 * Only valid if \c type is \ref TU_DEC_TWO_SUM or \ref TU_DEC_THREE_SUM.
 */
TU_EXPORT
int TUgetDecRankLowerLeft(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the rank of the top-right submatrix of this node of the tree.
 *
 * Only valid if \c type is \ref TU_DEC_TWO_SUM or \ref TU_DEC_THREE_SUM.
 */
TU_EXPORT
int TUgetDecRankTopRight(
  TU_DEC* dec /**< Decomposition tree */
);

bool TUtestBinaryRegularSequentiallyConnected(TU* tu, TU_CHAR_MATRIX* matrix,
  TU_CHAR_MATRIX* transpose, TU_DEC** decomposition, bool notGraphic, bool notCographic)
{
  return true;
}

int compareComponents(const void* a, const void* b)
{
  return ((TU_ONESUM_COMPONENT*)a)->matrix->numNonzeros -
    ((TU_ONESUM_COMPONENT*)b)->matrix->numNonzeros;
}

bool TUtestBinaryRegularLabeled(TU* tu, TU_CHAR_MATRIX* matrix, int* rowLabels, int* columnLabels,
  TU_DEC** decomposition)
{
  bool isRegular = true;
  int numComponents;
  TU_ONESUM_COMPONENT* components;

  assert(tu);
  assert(matrix);

  /* Check entries. */

  assert(TUisBinaryChar(tu, matrix, NULL));

  /* Perform 1-sum decomposition. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  /* Sort components by number of nonzeros. */
  TU_ONESUM_COMPONENT** orderedComponents = NULL;
  TUallocStackArray(tu, &orderedComponents, numComponents);
  for (int comp = 0; comp < numComponents; ++comp)
    orderedComponents[comp] = &components[comp];
  qsort(orderedComponents, numComponents, sizeof(TU_ONESUM_COMPONENT*), &compareComponents);

  /* Test regularity for each component. */

  if (numComponents > 0 && decomposition)
  {
    TUallocBlock(tu, decomposition);
    (*decomposition)->numChildren = numComponents;
    TUallocBlockArray(tu, &(*decomposition)->children, numComponents);
    (*decomposition)->graph = NULL;
    (*decomposition)->cograph = NULL;
    TUallocBlockArray(tu, &(*decomposition)->rowLabels, matrix->numRows);
    TUallocBlockArray(tu, &(*decomposition)->columnLabels, matrix->numColumns);
    for (int row = 0; row < (*decomposition)->matrix->numRows; ++row)
      (*decomposition)->rowLabels[row] = rowLabels[row];
    for (int column = 0; column < (*decomposition)->matrix->numColumns; ++column)
      (*decomposition)->columnLabels[column] = columnLabels[column];
    (*decomposition)->flags = TU_DEC_ONE_SUM;

    TUcopyCharMatrix(tu, matrix, &(*decomposition)->matrix);
    (*decomposition)->transpose = NULL;
  }

  for (int i = 0; i < numComponents; ++i)
  {
    int comp = (orderedComponents[i] - components) / sizeof(TU_ONESUM_COMPONENT*);
    assert(comp >= 0);
    assert(comp < numComponents);

    TU_DEC* compDecomposition = NULL;

    bool compRegular = TUtestBinaryRegularSequentiallyConnected(tu,
      (TU_CHAR_MATRIX*) &components[comp].matrix, (TU_CHAR_MATRIX*) &components[comp].transpose,
      decomposition ? &compDecomposition : NULL, false, false);

    isRegular = isRegular && compRegular;

    if (!decomposition && !compRegular)
      goto cleanup;

    if (decomposition)
    {
      assert(compRegular);
    }
  }

cleanup:
  TUfreeStackArray(tu, &orderedComponents);

  for (int comp = 0; comp < numComponents; ++comp)
  {
    if (decomposition)
    {
//       (*decomposition)->
    }
    else
    {
      TUfreeCharMatrix(tu, (TU_CHAR_MATRIX**) &components[comp].matrix);
      TUfreeCharMatrix(tu, (TU_CHAR_MATRIX**) &components[comp].transpose);
    }
    TUfreeBlockArray(tu, &components[comp].rowsToOriginal);
    TUfreeBlockArray(tu, &components[comp].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  return isRegular;
}


bool TUtestBinaryRegular(TU* tu, TU_CHAR_MATRIX* matrix, TU_DEC** decomposition)
{
  bool result;
  int* rowLabels = NULL;
  int* columnLabels = NULL;

  assert(tu);
  assert(matrix);

  TUallocStackArray(tu, &rowLabels, matrix->numRows);
  TUallocStackArray(tu, &columnLabels, matrix->numColumns);

  for (int row = 0; row < matrix->numRows; ++row)
    rowLabels[row] = -1 - row;
  for (int column = 0; column < matrix->numColumns; ++column)
    columnLabels[column] = 1 + column;

  result = TUtestBinaryRegularLabeled(tu, matrix, rowLabels, columnLabels, decomposition);

  TUfreeStackArray(tu, &columnLabels);
  TUfreeStackArray(tu, &rowLabels);

  return result;
}

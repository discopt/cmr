#include <tu/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "regular_internal.h"
#include "env_internal.h"
#include "one_sum.h"

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
  dec->flags = 0;
  dec->children = NULL;
  dec->numChildren = 0;
  dec->graph = NULL;
  dec->cograph = NULL;
}

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
int TUgetDecRankLowerLeft(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the rank of the top-right submatrix of this node of the tree.
 *
 * Only valid if \c type is \ref TU_DEC_TWO_SUM or \ref TU_DEC_THREE_SUM.
 */
int TUgetDecRankTopRight(
  TU_DEC* dec /**< Decomposition tree */
);


bool TUregularSequentiallyConnected(TU* tu, TU_DEC* decomposition, bool certify, bool notGraphic,
  bool notCographic)
{
  assert(tu);
  assert(decomposition);
  assert(decomposition->matrix);
  assert(decomposition->transpose);

  return true;
}

int compareOneSumComponents(const void* a, const void* b)
{
  return ((TU_ONESUM_COMPONENT*)a)->matrix->numNonzeros -
    ((TU_ONESUM_COMPONENT*)b)->matrix->numNonzeros;
}

int TUregularDecomposeOneSum(TU* tu, TU_CHAR_MATRIX* matrix, int* rowLabels, int* columnLabels,
  TU_DEC** pdecomposition)
{
  assert(tu);
  assert(matrix);
  assert(TUisBinaryChar(tu, matrix, NULL));
  assert(pdecomposition);

  TUcreateDec(tu, pdecomposition);
  TU_DEC* decomposition = *pdecomposition;

  /* Perform 1-sum decomposition. */

  int numComponents;
  TU_ONESUM_COMPONENT* components = NULL;
  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  if (numComponents <= 1)
  {
    decomposition->matrix = (TU_CHAR_MATRIX*) components[0].matrix;
    decomposition->transpose = (TU_CHAR_MATRIX*) components[0].transpose;
    if (rowLabels)
    {
      TUallocBlockArray(tu, &decomposition->rowLabels, matrix->numRows);
      for (int row = 0; row < matrix->numRows; ++row)
        decomposition->rowLabels[row] = rowLabels[components[0].rowsToOriginal[row]];
    }
    if (columnLabels)
    {
      TUallocBlockArray(tu, &decomposition->columnLabels, matrix->numColumns);
      for (int column = 0; column < matrix->numColumns; ++column)
        decomposition->columnLabels[column] = columnLabels[components[0].columnsToOriginal[column]];
    }
  }
  else
  {
    decomposition->flags = TU_DEC_ONE_SUM;

    /* Copy matrix and labels to node and compute transpose. */
    TUcopyCharMatrix(tu, matrix, &decomposition->matrix);
    TUtransposeCharMatrix(tu, matrix, &decomposition->transpose);
    if (rowLabels)
    {
      TUallocBlockArray(tu, &decomposition->rowLabels, matrix->numRows);
      for (int row = 0; row < matrix->numRows; ++row)
        decomposition->rowLabels[row] = rowLabels[row];
    }
    if (columnLabels)
    {
      TUallocBlockArray(tu, &decomposition->columnLabels, matrix->numColumns);
      for (int column = 0; column < matrix->numColumns; ++column)
        decomposition->columnLabels[column] = columnLabels[column];
    }

    /* Sort components by number of nonzeros. */
    TU_ONESUM_COMPONENT** orderedComponents = NULL;
    TUallocStackArray(tu, &orderedComponents, numComponents);
    for (int comp = 0; comp < numComponents; ++comp)
      orderedComponents[comp] = &components[comp];
    qsort(orderedComponents, numComponents, sizeof(TU_ONESUM_COMPONENT*), &compareOneSumComponents);

    /* Initialize child nodes */
    decomposition->numChildren = numComponents;
    TUallocBlockArray(tu, &decomposition->children, numComponents);

    for (int i = 0; i < numComponents; ++i)
    {
      int comp = (orderedComponents[i] - components) / sizeof(TU_ONESUM_COMPONENT*);
      TUcreateDec(tu, &decomposition->children[i]);
      TU_DEC* child = decomposition->children[i];
      child->matrix = (TU_CHAR_MATRIX*) components[comp].matrix;
      child->transpose = (TU_CHAR_MATRIX*) components[comp].transpose;
      if (rowLabels)
      {
        TUallocBlockArray(tu, &child->rowLabels, child->matrix->numRows);
        for (int row = 0; row < child->matrix->numRows; ++row)
          child->rowLabels[row] = rowLabels[components[comp].rowsToOriginal[row]];
      }
      if (columnLabels)
      {
        TUallocBlockArray(tu, &child->columnLabels, child->matrix->numColumns);
        for (int column = 0; column < child->matrix->numColumns; ++column)
          child->columnLabels[column] = columnLabels[components[comp].columnsToOriginal[column]];
      }
    }

    TUfreeStackArray(tu, &orderedComponents);
  }

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TUfreeBlockArray(tu, &components[comp].rowsToOriginal);
    TUfreeBlockArray(tu, &components[comp].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  return numComponents;
}

bool TUregularTest(TU* tu, TU_CHAR_MATRIX* matrix, int* rowLabels, int* columnLabels,
  TU_DEC** pdecomposition)
{
  bool isRegular = true;
  bool certify = pdecomposition != NULL;

  assert(tu);
  assert(matrix);

  /* Perform a 1-sum decomposition. */
  TU_DEC* decomposition = NULL;

  int numChildren = TUregularDecomposeOneSum(tu, matrix, rowLabels, columnLabels, &decomposition);
  if (certify)
    *pdecomposition = decomposition;
  if (numChildren <= 1)
  {
    isRegular = TUregularSequentiallyConnected(tu, decomposition, certify, false, false);
  }
  else
  {
    if (certify)
      decomposition->flags = TU_DEC_ONE_SUM | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC | TU_DEC_REGULAR;

    for (int child = 0; child < numChildren; ++child)
    {
      bool result = TUregularSequentiallyConnected(tu, decomposition->children[child], certify,
        false, false);
      isRegular = isRegular && result;
      if (certify)
      {
        decomposition->flags &= (decomposition->children[child]->flags
          & (TU_DEC_REGULAR | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC));
      }
      else if (!result)
        break;
    }
  }

  if (!pdecomposition)
    TUfreeDec(tu, &decomposition);

  return isRegular;
}

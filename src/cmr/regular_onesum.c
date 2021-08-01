#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "regular_internal.h"
#include "env_internal.h"
#include "one_sum.h"

int compareOneSumComponents(const void* a, const void* b)
{
  return ((TU_ONESUM_COMPONENT*)a)->matrix->numNonzeros -
    ((TU_ONESUM_COMPONENT*)b)->matrix->numNonzeros;
}

int TUregularDecomposeOneSum(TU* tu, TU_CHRMAT* matrix, int* rowLabels, int* columnLabels,
  TU_DEC** pdec, bool constructDecomposition)
{
  assert(tu);
  assert(matrix);
  assert(TUisTernaryChr(tu, matrix, NULL));
  assert(pdec);

  TUcreateDec(tu, pdec);
  TU_DEC* dec = *pdec;

  /* Perform 1-sum decomposition. */

  int numComponents;
  TU_ONESUM_COMPONENT* components = NULL;
  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  if (numComponents <= 1)
  {
    dec->matrix = (TU_CHRMAT*) components[0].matrix;
    if (constructDecomposition)
      dec->transpose = (TU_CHRMAT*) components[0].transpose;
    else
      TUchrmatFree(tu, (TU_CHRMAT**) &components[0].transpose);
    if (rowLabels)
    {
      TUallocBlockArray(tu, &dec->rowLabels, matrix->numRows);
      for (int row = 0; row < matrix->numRows; ++row)
        dec->rowLabels[row] = rowLabels[components[0].rowsToOriginal[row]];
    }
    if (columnLabels)
    {
      TUallocBlockArray(tu, &dec->columnLabels, matrix->numColumns);
      for (int column = 0; column < matrix->numColumns; ++column)
        dec->columnLabels[column] = columnLabels[components[0].columnsToOriginal[column]];
    }
  }
  else
  {
    dec->flags = TU_DEC_ONE_SUM;

    /* Copy matrix and labels to node and compute transpose. */
    if (constructDecomposition)
    {
      TUchrmatCopy(tu, matrix, &dec->matrix);
      TUchrmatTranspose(tu, matrix, &dec->transpose);
    }
    if (rowLabels)
    {
      TUallocBlockArray(tu, &dec->rowLabels, matrix->numRows);
      for (int row = 0; row < matrix->numRows; ++row)
        dec->rowLabels[row] = rowLabels[row];
    }
    if (columnLabels)
    {
      TUallocBlockArray(tu, &dec->columnLabels, matrix->numColumns);
      for (int column = 0; column < matrix->numColumns; ++column)
        dec->columnLabels[column] = columnLabels[column];
    }

    /* Sort components by number of nonzeros. */
    TU_ONESUM_COMPONENT** orderedComponents = NULL;
    TUallocStackArray(tu, &orderedComponents, numComponents);
    for (int comp = 0; comp < numComponents; ++comp)
      orderedComponents[comp] = &components[comp];
    qsort(orderedComponents, numComponents, sizeof(TU_ONESUM_COMPONENT*), &compareOneSumComponents);

    /* Initialize child nodes */
    dec->numChildren = numComponents;
    TUallocBlockArray(tu, &dec->children, numComponents);

    for (int i = 0; i < numComponents; ++i)
    {
      int comp = (orderedComponents[i] - components) / sizeof(TU_ONESUM_COMPONENT*);
      TUcreateDec(tu, &dec->children[i]);
      TU_DEC* child = dec->children[i];
      child->matrix = (TU_CHRMAT*) components[comp].matrix;
      if (constructDecomposition)
        child->transpose = (TU_CHRMAT*) components[comp].transpose;
      else
        TUchrmatFree(tu, (TU_CHRMAT**) &components[comp].transpose);
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

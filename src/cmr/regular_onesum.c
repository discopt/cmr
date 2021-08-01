#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "regular_internal.h"
#include "env_internal.h"
#include "one_sum.h"

int compareOneSumComponents(const void* a, const void* b)
{
  return ((CMR_ONESUM_COMPONENT*)a)->matrix->numNonzeros -
    ((CMR_ONESUM_COMPONENT*)b)->matrix->numNonzeros;
}

int CMRregularDecomposeOneSum(CMR* cmr, CMR_CHRMAT* matrix, int* rowLabels, int* columnLabels,
  CMR_TU_DEC** pdec, bool constructDecomposition)
{
  assert(cmr);
  assert(matrix);
  assert(CMRisTernaryChr(cmr, matrix, NULL));
  assert(pdec);

  CMRcreateDec(cmr, pdec);
  CMR_TU_DEC* dec = *pdec;

  /* Perform 1-sum decomposition. */

  int numComponents;
  CMR_ONESUM_COMPONENT* components = NULL;
  decomposeOneSum(cmr, (CMR_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  if (numComponents <= 1)
  {
    dec->matrix = (CMR_CHRMAT*) components[0].matrix;
    if (constructDecomposition)
      dec->transpose = (CMR_CHRMAT*) components[0].transpose;
    else
      CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[0].transpose);
    if (rowLabels)
    {
      CMRallocBlockArray(cmr, &dec->rowLabels, matrix->numRows);
      for (int row = 0; row < matrix->numRows; ++row)
        dec->rowLabels[row] = rowLabels[components[0].rowsToOriginal[row]];
    }
    if (columnLabels)
    {
      CMRallocBlockArray(cmr, &dec->columnLabels, matrix->numColumns);
      for (int column = 0; column < matrix->numColumns; ++column)
        dec->columnLabels[column] = columnLabels[components[0].columnsToOriginal[column]];
    }
  }
  else
  {
    dec->flags = CMR_TU_DEC_ONE_SUM;

    /* Copy matrix and labels to node and compute transpose. */
    if (constructDecomposition)
    {
      CMRchrmatCopy(cmr, matrix, &dec->matrix);
      CMRchrmatTranspose(cmr, matrix, &dec->transpose);
    }
    if (rowLabels)
    {
      CMRallocBlockArray(cmr, &dec->rowLabels, matrix->numRows);
      for (int row = 0; row < matrix->numRows; ++row)
        dec->rowLabels[row] = rowLabels[row];
    }
    if (columnLabels)
    {
      CMRallocBlockArray(cmr, &dec->columnLabels, matrix->numColumns);
      for (int column = 0; column < matrix->numColumns; ++column)
        dec->columnLabels[column] = columnLabels[column];
    }

    /* Sort components by number of nonzeros. */
    CMR_ONESUM_COMPONENT** orderedComponents = NULL;
    CMRallocStackArray(cmr, &orderedComponents, numComponents);
    for (int comp = 0; comp < numComponents; ++comp)
      orderedComponents[comp] = &components[comp];
    qsort(orderedComponents, numComponents, sizeof(CMR_ONESUM_COMPONENT*), &compareOneSumComponents);

    /* Initialize child nodes */
    dec->numChildren = numComponents;
    CMRallocBlockArray(cmr, &dec->children, numComponents);

    for (int i = 0; i < numComponents; ++i)
    {
      int comp = (orderedComponents[i] - components) / sizeof(CMR_ONESUM_COMPONENT*);
      CMRcreateDec(cmr, &dec->children[i]);
      CMR_TU_DEC* child = dec->children[i];
      child->matrix = (CMR_CHRMAT*) components[comp].matrix;
      if (constructDecomposition)
        child->transpose = (CMR_CHRMAT*) components[comp].transpose;
      else
        CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[comp].transpose);
      if (rowLabels)
      {
        CMRallocBlockArray(cmr, &child->rowLabels, child->matrix->numRows);
        for (int row = 0; row < child->matrix->numRows; ++row)
          child->rowLabels[row] = rowLabels[components[comp].rowsToOriginal[row]];
      }
      if (columnLabels)
      {
        CMRallocBlockArray(cmr, &child->columnLabels, child->matrix->numColumns);
        for (int column = 0; column < child->matrix->numColumns; ++column)
          child->columnLabels[column] = columnLabels[components[comp].columnsToOriginal[column]];
      }
    }

    CMRfreeStackArray(cmr, &orderedComponents);
  }

  for (int comp = 0; comp < numComponents; ++comp)
  {
    CMRfreeBlockArray(cmr, &components[comp].rowsToOriginal);
    CMRfreeBlockArray(cmr, &components[comp].columnsToOriginal);
  }
  CMRfreeBlockArray(cmr, &components);

  return numComponents;
}

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "regular_internal.h"
#include "env_internal.h"

CMR_TU_DEC* CMRregularDecomposeSimpleSums(CMR* cmr, CMR_TU_DEC* dec, bool unitVectors, bool paths, bool constructDecomposition)
{
  assert(cmr);
  assert(dec);
  assert(CMRisTernaryChr(cmr, dec->matrix, NULL));

  CMR_TU_DEC* result = dec;
  CMR_CHRMAT* matrix = dec->matrix;
  CMR_CHRMAT* transpose = dec->transpose;

  int* rowNonzeros = NULL;
  int* columnNonzeros = NULL;
  CMRallocStackArray(cmr, &rowNonzeros, matrix->numRows);
  CMRallocStackArray(cmr, &columnNonzeros, matrix->numColumns);

  /* Ensure that we have the transpose available. */

  if (constructDecomposition && result->transpose == NULL)
    CMRchrmatTranspose(cmr, result->matrix, &result->transpose);

  if (unitVectors)
  {
    int* queue = NULL;
    int queueBegin = 0;
    int queueEnd = 0;
    CMRallocStackArray(cmr, &queue, matrix->numRows + matrix->numColumns);

    for (int row = 0; row < matrix->numRows; ++row)
    {
      rowNonzeros[row] =
        (row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros)
        - matrix->rowStarts[row];
      if (rowNonzeros[row] == 1)
      {
        queue[queueEnd] = -1 - row;
        queueEnd++;
      }
    }

    for (int column = 0; column < matrix->numColumns; ++column)
    {
      columnNonzeros[column] =
        (column + 1 < matrix->numColumns ? transpose->rowStarts[column + 1] : matrix->numNonzeros)
        - transpose->rowStarts[column];
      if (columnNonzeros[column] == 1)
      {
        queue[queueEnd] = column;
        queueEnd++;
      }
    }

    while (queueBegin < queueEnd)
    {
      int element = queue[queueBegin];
      result->numChildren = 2;
      result->flags = constructDecomposition ? CMR_TU_DEC_TWO_SUM : CMR_TU_DEC_MULTI_SUM;
      CMRallocBlockArray(cmr, &result->children, 2);
      CMRcreateDec(cmr, &result->children[0]);
      CMRcreateDec(cmr, &result->children[1]);

      /* Create children: First is 1x2 or 2x1 matrix, second is without the unit vector. */

      result->children[0]->flags = CMR_TU_DEC_PROCESSED | CMR_TU_DEC_GRAPHIC | CMR_TU_DEC_COGRAPHIC
        | CMR_TU_DEC_REGULAR;
      if (constructDecomposition)
      {
        if (element < 0)
        {
          int row = -1 - element;
          CMRchrmatCreate(cmr, &result->children[0]->matrix, 2, 1, 2);
          result->children[1]->matrix->rowStarts[0] = 0;
          result->children[1]->matrix->rowStarts[1] = 1;
          result->children[1]->matrix->rowStarts[2] = 2;
          result->children[1]->matrix->entryColumns[0] = 0;
          result->children[1]->matrix->entryColumns[1] = 0;
          result->children[1]->matrix->entryValues[0] = matrix->entryValues[matrix->rowStarts[row]];
          result->children[1]->matrix->entryValues[1] = 1;
          int column = matrix->entryColumns[matrix->rowStarts[row]];
          int begin = transpose->rowStarts[column];
          int end = column < transpose->numRows ? transpose->rowStarts[column + 1]
            : transpose->numNonzeros;
          int otherRow = -1;
          for (int entry = begin; entry < end; ++entry)
          {
            if (transpose->entryColumns[entry] != row)
            {
              result->children[1]->matrix->entryValues[1] = transpose->entryValues[entry];
              otherRow = transpose->entryColumns[entry];
              break;
            }
          }
          assert(otherRow >= 0);
          rowNonzeros[row] = 0;
          columnNonzeros[column]--;
          assert(columnNonzeros[column] >= 1);
          if (columnNonzeros[column] == 1)
          {
            queue[queueEnd] = column;
            queueEnd++;
          }
          if (dec->rowLabels && dec->columnLabels)
          {
            CMRallocBlockArray(cmr, &result->children[0]->rowLabels, 2);
            CMRallocBlockArray(cmr, &result->children[0]->columnLabels, 1);
            result->children[0]->rowLabels[0] = dec->rowLabels[row];
            result->children[0]->rowLabels[1] = dec->rowLabels[otherRow];
            result->children[0]->columnLabels[0] = dec->columnLabels[column];
          }

          /* Create matrix of 2nd child. */

          CMRchrmatCreate(cmr, &result->children[1]->matrix, result->matrix->numRows - 1,
            result->matrix->numColumns, result->matrix->numNonzeros - 1);
          int entry = 0;
          for (int r = 0; r < result->matrix->numRows; ++r)
          {
            result->children[1]->matrix->rowStarts[r] = entry;
            int begin, end, offset;
            if (r < row)
            {
              begin = result->matrix->rowStarts[r];
              end = r + 1 < result->matrix->numRows ? result->matrix->rowStarts[r + 1]
                : result->matrix->numNonzeros;
              offset = 0;
            }
            else if (r > row)
            {
              begin = result->matrix->rowStarts[r + 1];
              end = r + 2 < result->matrix->numRows ? result->matrix->rowStarts[r + 2]
                : result->matrix->numNonzeros;
              offset = 1;
            }
            else
              continue;
            for (int e = begin; e < end; ++e)
            {
              result->children[1]->matrix->entryColumns[e] = result->matrix->entryColumns[e + offset];
              result->children[1]->matrix->entryValues[e] = result->matrix->entryValues[e + offset];
            }
          }
          CMRchrmatTranspose(cmr, result->children[1]->matrix, &result->children[1]->transpose);
        }
        else
        {
          int column = element;
          CMRchrmatCreate(cmr, &result->children[0]->matrix, 1, 2, 2);
          result->children[1]->matrix->rowStarts[0] = 0;
          result->children[1]->matrix->rowStarts[1] = 2;
          result->children[1]->matrix->entryColumns[0] = 0;
          result->children[1]->matrix->entryColumns[1] = 1;
          result->children[1]->matrix->entryValues[0] = transpose->entryValues[transpose->rowStarts[column]];
          result->children[1]->matrix->entryValues[1] = 1;
          int row = transpose->entryColumns[transpose->rowStarts[column]];
          int begin = matrix->rowStarts[row];
          int end = row < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;
          int otherColumn = -1;
          for (int entry = begin; entry < end; ++entry)
          {
            if (matrix->entryColumns[entry] != column)
            {
              result->children[1]->matrix->entryValues[1] = matrix->entryValues[entry];
              otherColumn = matrix->entryColumns[entry];
              break;
            }
          }
          assert(otherColumn >= 0);
          columnNonzeros[column] = 0;
          rowNonzeros[row]--;
          assert(rowNonzeros[row] >= 1);
          if (rowNonzeros[row] == 1)
          {
            queue[queueEnd] = -1 - row;
            queueEnd++;
          }
          if (dec->rowLabels && dec->columnLabels)
          {
            CMRallocBlockArray(cmr, &result->children[0]->rowLabels, 1);
            CMRallocBlockArray(cmr, &result->children[0]->columnLabels, 2);
            result->children[0]->columnLabels[0] = dec->columnLabels[column];
            result->children[0]->columnLabels[1] = dec->columnLabels[otherColumn];
            result->children[0]->rowLabels[0] = dec->rowLabels[row];
          }

          /* Create matrix of 2nd child. */

          CMRchrmatCreate(cmr, &result->children[1]->transpose, result->transpose->numRows - 1,
            result->transpose->numColumns, result->transpose->numNonzeros - 1);
          int entry = 0;
          for (int r = 0; r < result->transpose->numRows; ++r)
          {
            result->children[1]->transpose->rowStarts[r] = entry;
            int begin, end, offset;
            if (r < row)
            {
              begin = result->transpose->rowStarts[r];
              end = r + 1 < result->transpose->numRows ? result->transpose->rowStarts[r + 1]
                : result->transpose->numNonzeros;
              offset = 0;
            }
            else if (r > row)
            {
              begin = result->transpose->rowStarts[r + 1];
              end = r + 2 < result->transpose->numRows ? result->transpose->rowStarts[r + 2]
                : result->transpose->numNonzeros;
              offset = 1;
            }
            else
              continue;
            for (int e = begin; e < end; ++e)
            {
              result->children[1]->transpose->entryColumns[e] = result->transpose->entryColumns[e + offset];
              result->children[1]->transpose->entryValues[e] = result->transpose->entryValues[e + offset];
            }
          }
          CMRchrmatTranspose(cmr, result->children[1]->transpose, &result->children[1]->matrix);
        }
        CMRchrmatTranspose(cmr, result->children[0]->matrix, &result->children[0]->transpose);
      }

      /* Make second child the current node. */

      result = result->children[1];
    }

    CMRfreeStackArray(cmr, &queue);
  }

  CMRfreeStackArray(cmr, &columnNonzeros);
  CMRfreeStackArray(cmr, &rowNonzeros);

  return result;
}

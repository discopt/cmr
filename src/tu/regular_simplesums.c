#include <tu/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "regular_internal.h"
#include "env_internal.h"

TU_DEC* TUregularDecomposeSimpleSums(TU* tu, TU_DEC* decomposition, bool unitVectors, bool paths,
  bool constructDecomposition)
{
  assert(tu);
  assert(decomposition);
  assert(TUisTernaryChr(tu, decomposition->matrix, NULL));

  TU_DEC* result = decomposition;
  TU_CHRMAT* matrix = decomposition->matrix;
  TU_CHRMAT* transpose = decomposition->transpose;

  int* rowNonzeros = NULL;
  int* columnNonzeros = NULL;
  TUallocStackArray(tu, &rowNonzeros, matrix->numRows);
  TUallocStackArray(tu, &columnNonzeros, matrix->numColumns);

  /* Ensure that we have the transpose available. */

  if (constructDecomposition && result->transpose == NULL)
    TUchrmatTranspose(tu, result->matrix, &result->transpose);

  if (unitVectors)
  {
    int* queue = NULL;
    int queueBegin = 0;
    int queueEnd = 0;
    TUallocStackArray(tu, &queue, matrix->numRows + matrix->numColumns);

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
      result->flags = constructDecomposition ? TU_DEC_TWO_SUM : TU_DEC_MULTI_SUM;
      TUallocBlockArray(tu, &result->children, 2);
      TUcreateDec(tu, &result->children[0]);
      TUcreateDec(tu, &result->children[1]);

      /* Create children: First is 1x2 or 2x1 matrix, second is without the unit vector. */

      result->children[0]->flags = TU_DEC_PROCESSED | TU_DEC_GRAPHIC | TU_DEC_COGRAPHIC
        | TU_DEC_REGULAR;
      if (constructDecomposition)
      {
        if (element < 0)
        {
          int row = -1 - element;
          TUchrmatCreate(tu, &result->children[0]->matrix, 2, 1, 2);
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
          if (decomposition->rowLabels && decomposition->columnLabels)
          {
            TUallocBlockArray(tu, &result->children[0]->rowLabels, 2);
            TUallocBlockArray(tu, &result->children[0]->columnLabels, 1);
            result->children[0]->rowLabels[0] = decomposition->rowLabels[row];
            result->children[0]->rowLabels[1] = decomposition->rowLabels[otherRow];
            result->children[0]->columnLabels[0] = decomposition->columnLabels[column];
          }

          /* Create matrix of 2nd child. */

          TUchrmatCreate(tu, &result->children[1]->matrix, result->matrix->numRows - 1,
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
          TUchrmatTranspose(tu, result->children[1]->matrix, &result->children[1]->transpose);
        }
        else
        {
          int column = element;
          TUchrmatCreate(tu, &result->children[0]->matrix, 1, 2, 2);
          result->children[1]->matrix->rowStarts[0] = 0;
          result->children[1]->matrix->rowStarts[1] = 2;
          result->children[1]->matrix->entryColumns[0] = 0;
          result->children[1]->matrix->entryColumns[1] = 1;
          result->children[1]->matrix->entryValues[0] 
            = transpose->entryValues[transpose->rowStarts[column]];
          result->children[1]->matrix->entryValues[1] = 1;
          int row = transpose->entryColumns[transpose->rowStarts[column]];
          int begin = matrix->rowStarts[row];
          int end = row < matrix->numRows ? matrix->rowStarts[row + 1] 
            : matrix->numNonzeros;
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
          if (decomposition->rowLabels && decomposition->columnLabels)
          {
            TUallocBlockArray(tu, &result->children[0]->rowLabels, 1);
            TUallocBlockArray(tu, &result->children[0]->columnLabels, 2);
            result->children[0]->columnLabels[0] = decomposition->columnLabels[column];
            result->children[0]->columnLabels[1] = decomposition->columnLabels[otherColumn];
            result->children[0]->rowLabels[0] = decomposition->rowLabels[row];
          }

          /* Create matrix of 2nd child. */

          TUchrmatCreate(tu, &result->children[1]->transpose, result->transpose->numRows - 1,
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
          TUchrmatTranspose(tu, result->children[1]->transpose, &result->children[1]->matrix);
        }
        TUchrmatTranspose(tu, result->children[0]->matrix, &result->children[0]->transpose);
      }

      /* Make second child the current node. */
      
      result = result->children[1];
    }

    TUfreeStackArray(tu, &queue);
  }

  TUfreeStackArray(tu, &columnNonzeros);
  TUfreeStackArray(tu, &rowNonzeros);

  return result;
}

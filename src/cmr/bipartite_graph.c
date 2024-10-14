// #define CMR_DEBUG /** Uncomment to debug this file. */

#include "bipartite_graph.h"

#include "env_internal.h"

CMR_ERROR CMRchrmatSubmatrixBipartitePath(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose, int* rowsGroup,
  int* columnsGroup, bool* pconnected, CMR_ELEMENT* ppathSource, CMR_ELEMENT* ppathTarget, CMR_ELEMENT* rowsPredecessor,
  CMR_ELEMENT* columnsPredecessor, int* psum)
{
  assert(cmr);
  assert(matrix);
  assert(transpose);
  assert(rowsGroup);
  assert(columnsGroup);

  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;

  CMRdbgMsg(12, "CMRchrmatSubmatrixBipartitePath() called.\n");

  bool localRowsPredecessor = rowsPredecessor == NULL;
  if (localRowsPredecessor)
    CMR_CALL( CMRallocStackArray(cmr, &rowsPredecessor, numRows) );
  bool localColumnsPredecessor = columnsPredecessor == NULL;
  if (localColumnsPredecessor)
    CMR_CALL( CMRallocStackArray(cmr, &columnsPredecessor, numColumns) );

  int8_t* rowsStatus = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowsStatus, numRows) );
  int8_t* columnsStatus = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnsStatus, numColumns) );
  int8_t* rowsPredecessorEdgeValue = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowsPredecessorEdgeValue, numRows) );
  int8_t* columnsPredecessorEdgeValue = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnsPredecessorEdgeValue, numColumns) );
  CMR_ELEMENT* queue = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &queue, numRows + numColumns) );

  bool connected = false;
  CMR_ELEMENT pathTarget = 0;
  for (int sourceGroup = 1; !connected; ++sourceGroup)
  {
    size_t firstQueued = 0;
    size_t beyondQueued = 0;
    for (size_t row = 0; row < numRows; ++row)
    {
      if (rowsGroup[row] == sourceGroup)
      {
        queue[beyondQueued++] = CMRrowToElement(row);
        rowsStatus[row] = 1;
        rowsPredecessor[row] = 0;
      }
      else if (rowsGroup[row] < 0)
        rowsStatus[row] = -1;
      else
        rowsStatus[row] = 0;
    }
    for (size_t column = 0; column < numColumns; ++column)
    {
      if (columnsGroup[column] == sourceGroup)
      {
        queue[beyondQueued++] = CMRcolumnToElement(column);
        columnsStatus[column] = 1;
        columnsPredecessor[column] = 0;
      }
      else if (columnsGroup[column] < 0)
        columnsStatus[column] = -1;
      else
        columnsStatus[column] = 0;
    }

    if (beyondQueued == 0)
      break;

    /** Run BFS. */
    while (firstQueued < beyondQueued)
    {
      CMR_ELEMENT current = queue[firstQueued++];

      CMRdbgMsg(12, "Processing %s.\n", CMRelementString(current, NULL));

      if (CMRelementIsRow(current))
      {
        size_t row = CMRelementToRowIndex(current);
        rowsStatus[row] = 2;

        size_t first = matrix->rowSlice[row];
        size_t beyond = matrix->rowSlice[row + 1];
        for (size_t e = first; e < beyond; ++e)
        {
          size_t column = matrix->entryColumns[e];
          if (columnsStatus[column])
            continue;

          columnsPredecessor[column] = current;
          columnsPredecessorEdgeValue[column] = matrix->entryValues[e];

          if (columnsGroup[column] > sourceGroup)
          {
            pathTarget = CMRcolumnToElement(column);
            firstQueued = beyondQueued;
            connected = true;
            break;
          }

          columnsStatus[column] = 1;
          queue[beyondQueued++] = CMRcolumnToElement(column);
        }
      }
      else
      {
        size_t column = CMRelementToColumnIndex(current);
        columnsStatus[column] = 2;

        size_t first = transpose->rowSlice[column];
        size_t beyond = transpose->rowSlice[column + 1];
        for (size_t e = first; e < beyond; ++e)
        {
          size_t row = transpose->entryColumns[e];
          if (rowsStatus[row])
            continue;

          rowsPredecessor[row] = current;
          rowsPredecessorEdgeValue[row] = transpose->entryValues[e];

          if (rowsGroup[row] > sourceGroup)
          {
            pathTarget = CMRrowToElement(row);
            firstQueued = beyondQueued;
            connected = true;
            break;
          }

          rowsStatus[row] = 1;
          queue[beyondQueued++] = CMRrowToElement(row);
        }
      }
    }
  }

  if (pconnected)
    *pconnected = connected;
  if (ppathTarget)
    *ppathTarget = pathTarget;

  CMR_ELEMENT current = pathTarget;
  int sum = 0;
  while (CMRelementIsValid(current))
  {
    if (CMRelementIsRow(current))
    {
      size_t row = CMRelementToRowIndex(current);
      CMR_ELEMENT pred = rowsPredecessor[row];
      if (!CMRelementIsValid(pred))
        break;
      sum += rowsPredecessorEdgeValue[row];
      current = pred;
    }
    else
    {
      size_t column = CMRelementToColumnIndex(current);
      CMR_ELEMENT pred = columnsPredecessor[column];
      if (!CMRelementIsValid(pred))
        break;
      sum += columnsPredecessorEdgeValue[column];
      current = pred;
    }
  }

  if (ppathSource)
    *ppathSource = current;
  if (psum)
    *psum = sum;

  /* Cleanup */

  CMR_CALL( CMRfreeStackArray(cmr, &queue) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnsPredecessorEdgeValue) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowsPredecessorEdgeValue) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnsStatus) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowsStatus) );
  if (localColumnsPredecessor)
    CMR_CALL( CMRfreeStackArray(cmr, &columnsPredecessor) );
  if (localRowsPredecessor)
    CMR_CALL( CMRfreeStackArray(cmr, &rowsPredecessor) );

  return CMR_OKAY;
}



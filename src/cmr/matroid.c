#define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/matroid.h>


#include "env_internal.h"
#include "matroid_internal.h"
#include "matrix_internal.h"
#include "listmatrix.h"

#include <assert.h>
#include <string.h>

static
CMR_ERROR computePivots(CMR* cmr, CMR_CHRMAT* matrix, size_t numPivots, size_t* pivotRows, size_t* pivotColumns,
  int8_t characteristic, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);
  assert(!numPivots || pivotRows);
  assert(!numPivots || pivotColumns);

#if defined(CMR_DEBUG)
  CMRdbgMsg(2, "Applying %zu pivots to %zux%zu matrix.\n", numPivots, matrix->numRows, matrix->numColumns);
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
#endif /* CMR_DEBUG */

  ListMat8* listmat = NULL;
  CMR_ERROR error = CMR_OKAY;
  size_t memNonzeros = 2 * matrix->numNonzeros + 256;
  CMR_CALL( CMRlistmat8Alloc(cmr, matrix->numRows, matrix->numColumns, memNonzeros, &listmat) );
  CMR_CALL( CMRlistmat8InitializeFromChrMatrix(cmr, listmat, matrix) );

  size_t* affectedRows = NULL;
  size_t numAffectedRows;
  CMR_CALL( CMRallocStackArray(cmr, &affectedRows, matrix->numRows) );
  int* denseColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseColumn, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    denseColumn[row] = INT_MIN;

  size_t* affectedColumns = NULL;
  size_t numAffectedColumns;
  CMR_CALL( CMRallocStackArray(cmr, &affectedColumns, matrix->numColumns) );
  int* denseRow = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseRow, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    denseRow[column] = INT_MIN;

  for (size_t pivot = 0; pivot < numPivots; ++pivot)
  {
    size_t pivotRow = pivotRows[pivot];
    size_t pivotColumn = pivotColumns[pivot];

    /* Compute the pivot value and all columns that are affected. */
    numAffectedColumns = 0;
    ListMat8Nonzero* head = &listmat->rowElements[pivotRow].head;
    for (ListMat8Nonzero* nz = head->right; nz != head; nz = nz->right)
    {
      if (denseRow[nz->column] == INT_MIN)
      {
        affectedColumns[numAffectedColumns++] = nz->column;
        denseRow[nz->column] = 0;
      }
      denseRow[nz->column] += nz->value;
    }

    int pivotValue = denseRow[pivotColumn];
    if (pivotValue == INT_MIN || (pivotValue % characteristic) == 0)
    {
      error = CMR_ERROR_INPUT;
      goto cleanup;
    }
    pivotValue %= characteristic;
    if (pivotValue < 0)
      pivotValue += characteristic;

    /* Compute all rows that are affected. */
    head = &listmat->columnElements[pivotColumn].head;
    numAffectedRows = 0;
    for (ListMat8Nonzero* nz = head->below; nz != head; nz = nz->below)
    {
      if (denseColumn[nz->row] == INT_MIN)
      {
        affectedRows[numAffectedRows++] = nz->row;
        denseColumn[nz->row] = 0;
      }
      denseColumn[nz->row] += nz->value;
    }

    /* Apply the pivot to all other entries first. */
    for (size_t r = 0; r < numAffectedRows; ++r)
    {
      size_t row = affectedRows[r];
      if (row == pivotRow)
        continue;

      int rowValue = denseColumn[row] % characteristic;
      if (rowValue == 0)
        continue;
      if (rowValue < 0)
        rowValue += characteristic;

      for (size_t c = 0; c < numAffectedColumns; ++c)
      {
        size_t column = affectedColumns[c];
        if (column == pivotColumn)
          continue;

        int columnValue = denseRow[column] % characteristic;
        if (columnValue == 0)
          continue;
        if (columnValue < 0)
          columnValue += characteristic;

        CMR_CALL( CMRlistmat8Insert(cmr, listmat, row, column, -pivotValue * rowValue * columnValue, 0, NULL) );
      }
    }

#if defined(CMR_DEBUG)
    CMRdbgMsg(4, "Pivot at r%zu,c%zu; after updating non-pivot rows/columns:\n", pivotRow+1, pivotColumn+1);
    CMRlistmat8PrintDense(cmr, listmat, stdout);
#endif /* CMR_DEBUG */

    /* Apply the pivot to the pivot row and clear the dense vector. */
    head = &listmat->rowElements[pivotRow].head;
    for (ListMat8Nonzero* nz = head->right; nz != head; nz = nz->right)
    {
      if (pivotValue == 2 && nz->column != pivotColumn)
        nz->value *= 2;
      denseRow[nz->column] = INT_MIN;
    }

    /* Apply the pivot to the pivot column and clear the dense vector. */
    head = &listmat->columnElements[pivotColumn].head;
    for (ListMat8Nonzero* nz = head->below; nz != head; nz = nz->below)
    {
      if (pivotValue == 2 && nz->row != pivotRow)
        nz->value *= 2;
      else if (nz->row == pivotRow)
      {
        nz->value *= -1;
      }
      denseColumn[nz->row] = INT_MIN;
    }

#if defined(CMR_DEBUG)
    CMRdbgMsg(4, "Pivot at r%zu,c%zu; after updating all rows/columns:\n", pivotRow+1, pivotColumn+1);
    CMRlistmat8PrintDense(cmr, listmat, stdout);
#endif /* CMR_DEBUG */

  }

  /* Extract all nonzeros and merge duplicates. */
  CMR_CALL( CMRchrmatCreate(cmr, presult, matrix->numRows, matrix->numColumns, memNonzeros) );
  CMR_CHRMAT* result = *presult;
  size_t entry = 0;
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    /* Compute all columns with nonzeros in the row. */
    numAffectedColumns = 0;
    ListMat8Nonzero* head = &listmat->rowElements[row].head;
    for (ListMat8Nonzero* nz = head->right; nz != head; nz = nz->right)
    {
      CMRdbgMsg(6, "Considering nonzero %d at r%zu,c%zu.\n", nz->value, nz->row+1, nz->column+1);
      if (denseRow[nz->column] == INT_MIN)
      {
        affectedColumns[numAffectedColumns++] = nz->column;
        denseRow[nz->column] = 0;
      }
      denseRow[nz->column] += nz->value;
    }

    if (entry + numAffectedColumns > result->numNonzeros)
    {
      size_t newNumNonzeros = 2 * result->numNonzeros + 256;
      CMR_CALL( CMRchrmatChangeNumNonzeros(cmr, result, newNumNonzeros) );
    }

    CMRdbgMsg(8, "#entries in row r%zu is <= %zu.\n", row+1, numAffectedColumns);

    result->rowSlice[row] = entry;
    for (size_t i = 0; i < numAffectedColumns; ++i)
    {
      size_t column = affectedColumns[i];
      int value = denseRow[column] % characteristic;
      CMRdbgMsg(6, "Aggregated value in r%zu,c%zu is %d.\n", row+1, column+1, value);
      denseRow[column] = INT_MIN;
      if (value == 0)
        continue;

      value = (value + characteristic) % characteristic;
      if (value == 2)
        value = -1;

      result->entryColumns[entry] = column;
      result->entryValues[entry] = value;
      ++entry;
    }
  }
  result->rowSlice[matrix->numRows] = entry;
  result->numNonzeros = entry;

  CMR_CALL( CMRchrmatSortNonzeros(cmr, result) );

cleanup:

  CMR_CALL( CMRfreeStackArray(cmr, &denseRow) );
  CMR_CALL( CMRfreeStackArray(cmr, &affectedColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &denseColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &affectedRows) );
  CMR_CALL( CMRlistmat8Free(cmr, &listmat) );

  return error;
}

CMR_ERROR CMRchrmatBinaryPivot(CMR* cmr, CMR_CHRMAT* matrix, size_t pivotRow, size_t pivotColumn, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, 1, &pivotRow, &pivotColumn, 2, presult) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatTernaryPivot(CMR* cmr, CMR_CHRMAT* matrix, size_t pivotRow, size_t pivotColumn, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, 1, &pivotRow, &pivotColumn, 3, presult) );

  return CMR_OKAY;
}



CMR_ERROR CMRchrmatBinaryPivots(CMR* cmr, CMR_CHRMAT* matrix, size_t numPivots, size_t* pivotRows, size_t* pivotColumns,
  CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(pivotRows);
  assert(pivotColumns);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, numPivots, pivotRows, pivotColumns, 2, presult) );

  return CMR_OKAY;
}

CMR_ERROR CMRchrmatTernaryPivots(CMR* cmr, CMR_CHRMAT* matrix, size_t numPivots, size_t* pivotRows,
  size_t* pivotColumns, CMR_CHRMAT** presult)
{
  assert(cmr);
  assert(matrix);
  assert(pivotRows);
  assert(pivotColumns);
  assert(presult);

  CMR_CALL( computePivots(cmr, matrix, numPivots, pivotRows, pivotColumns, 3, presult) );

  return CMR_OKAY;
}




CMR_ERROR CMRminorCreate(CMR* cmr, CMR_MINOR** pminor, size_t numPivots, CMR_SUBMAT* submatrix)
{
  assert(cmr);
  assert(pminor);
  assert(!*pminor);

  CMR_CALL( CMRallocBlock(cmr, pminor) );
  CMR_MINOR* minor = *pminor;
  minor->numPivots = numPivots;
  CMR_CALL( CMRallocBlockArray(cmr, &minor->pivotRows, numPivots) );
  CMR_CALL( CMRallocBlockArray(cmr, &minor->pivotColumns, numPivots) );
  minor->remainingSubmatrix = submatrix;

  return CMR_OKAY;
}

CMR_ERROR CMRminorFree(CMR* cmr, CMR_MINOR** pminor)
{
  assert(cmr);
  assert(pminor);

  if (!*pminor)
    return CMR_OKAY;

  CMR_MINOR* minor = *pminor;
  CMR_CALL( CMRfreeBlockArray(cmr, &minor->pivotRows) );
  CMR_CALL( CMRfreeBlockArray(cmr, &minor->pivotColumns) );
  CMR_CALL( CMRsubmatFree(cmr, &minor->remainingSubmatrix) );
  CMR_CALL( CMRfreeBlock(cmr, pminor) );

  return CMR_OKAY;
}

CMR_ERROR CMRminorWriteToStream(CMR* cmr, CMR_MINOR* minor, size_t numRows, size_t numColumns, FILE* stream)
{
  assert(cmr);
  assert(minor);
  assert(stream);

  CMR_CALL( CMRsubmatWriteToStream(cmr, minor->remainingSubmatrix, numRows, numColumns, stream) );

  return CMR_OKAY;
}


CMR_ERROR CMRminorWriteToFile(CMR* cmr, CMR_MINOR* minor, size_t numRows, size_t numColumns, const char* fileName)
{
  assert(cmr);
  assert(minor);

  FILE* stream;
  if (!fileName || !strcmp(fileName, "-"))
    stream = stdout;
  else
  {
    stream = fopen(fileName, "w");
    if (!stream)
      return CMR_ERROR_OUTPUT;
  }

  CMR_CALL( CMRminorWriteToStream(cmr, minor, numRows, numColumns, stream) );

  if (stream != stdout)
    fclose(stream);

  return CMR_OKAY;
}



bool CMRmatroiddecIsTernary(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->isTernary;
}

bool CMRmatroiddecThreeSumDistributedRanks(CMR_MATROID_DEC* dec)
{
  assert(dec);
  assert(dec->type == CMR_MATROID_DEC_TYPE_THREE_SUM);

  return dec->threesumFlags & CMR_MATROID_DEC_THREESUM_FLAG_DISTRIBUTED_RANKS;
}

bool CMRmatroiddecThreeSumConcentratedRank(CMR_MATROID_DEC* dec)
{
  assert(dec);
  assert(dec->type == CMR_MATROID_DEC_TYPE_THREE_SUM);

  return dec->threesumFlags & CMR_MATROID_DEC_THREESUM_FLAG_CONCENTRATED_RANK;
}

CMR_CHRMAT* CMRmatroiddecGetMatrix(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->matrix;
}

bool CMRmatroiddecHasTranspose(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->transpose != NULL;
}

CMR_CHRMAT* CMRmatroiddecGetTranspose(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->transpose;
}


size_t CMRmatroiddecNumChildren(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->numChildren;
}

CMR_MATROID_DEC* CMRmatroiddecChild(CMR_MATROID_DEC* dec, size_t childIndex)
{
  assert(dec);
  assert(childIndex < dec->numChildren);

  return dec->children[childIndex];
}

CMR_MATROID_DEC_TYPE CMRmatroiddecType(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->type;
}

int8_t CMRmatroiddecGraphicness(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->graphicness;
}

int8_t CMRmatroiddecCographicness(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->cographicness;
}

int8_t CMRmatroiddecRegularity(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->regularity;
}

size_t CMRmatroiddecNumRows(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->numRows;
}

size_t* CMRmatroiddecRowsParent(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->rowsParent;
}


size_t CMRmatroiddecNumColumns(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->numColumns;
}

CMR_ELEMENT * CMRmatroiddecColumnsRootElement(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->columnsRootElement;
}

size_t* CMRmatroiddecColumnsParent(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->columnsParent;
}

CMR_GRAPH* CMRmatroiddecGraph(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->graph;
}

CMR_GRAPH_EDGE* CMRmatroiddecGraphForest(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->graphForest;
}

size_t CMRmatroiddecGraphSizeForest(CMR_MATROID_DEC* dec)
{
  assert(dec);

  if (dec->matrix)
    return dec->numRows;
  else if (dec->transpose)
    return dec->numColumns;
  else if (dec->graph)
    return CMRgraphNumNodes(dec->graph) - 1;
  else
    return SIZE_MAX;
}

CMR_GRAPH_EDGE* CMRmatroiddecGraphCoforest(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->graphCoforest;
}

size_t CMRmatroiddecGraphSizeCoforest(CMR_MATROID_DEC* dec)
{
  assert(dec);

  if (dec->matrix)
    return dec->numColumns;
  else if (dec->transpose)
    return dec->numRows;
  else if (dec->graph)
    return CMRgraphNumEdges(dec->graph) + 1 - CMRgraphNumNodes(dec->graph);
  else
    return SIZE_MAX;
}

bool* CMRmatroiddecGraphArcsReversed(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->graphArcsReversed;
}

CMR_GRAPH* CMRmatroiddecCograph(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->cograph;
}

size_t CMRmatroiddecCographSizeForest(CMR_MATROID_DEC* dec)
{
  assert(dec);

  if (dec->matrix)
    return dec->numColumns;
  else if (dec->transpose)
    return dec->numRows;
  else if (dec->cograph)
    return CMRgraphNumNodes(dec->cograph) - 1;
  else
    return SIZE_MAX;
}

CMR_GRAPH_EDGE* CMRmatroiddecCographForest(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->cographForest;
}

size_t CMRmatroiddecCographSizeCoforest(CMR_MATROID_DEC* dec)
{
  assert(dec);

  if (dec->matrix)
    return dec->numRows;
  else if (dec->transpose)
    return dec->numColumns;
  else if (dec->cograph)
    return CMRgraphNumEdges(dec->cograph) + 1 - CMRgraphNumNodes(dec->cograph);
  else
    return SIZE_MAX;
}

CMR_GRAPH_EDGE* CMRmatroiddecCographCoforest(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->cographCoforest;
}

bool* CMRmatroiddecCographArcsReversed(CMR_MATROID_DEC* dec)
{
  assert(dec);

  return dec->cographArcsReversed;
}

CMR_ERROR CMRmatroiddecPrint(CMR* cmr, CMR_MATROID_DEC* dec, FILE* stream, size_t indent, bool printMatrices,
  bool printGraphs, bool printReductions)
{
  assert(cmr);
  assert(stream);

  /* Indent. */
  for (size_t i = 0; i < indent; ++i)
    fputc(' ', stream);

  if (!dec)
  {
    fprintf(stream, "<NULL>\n");
    return CMR_OKAY;
  }

  fprintf(stream, "%zux%zu ", dec->numRows, dec->numColumns);
  switch (dec->type)
  {
  case CMR_MATROID_DEC_TYPE_IRREGULAR:
    fprintf(stream, "irregular {");
  break;
  case CMR_MATROID_DEC_TYPE_UNKNOWN:
    fprintf(stream, "unknown {");
  break;
  case CMR_MATROID_DEC_TYPE_ONE_SUM:
    fprintf(stream, "1-sum with %zu children {", dec->numChildren);
  break;
  case CMR_MATROID_DEC_TYPE_TWO_SUM:
    fprintf(stream, "2-sum {");
    break;
  case CMR_MATROID_DEC_TYPE_THREE_SUM:
    fprintf(stream, "3-sum {");
  break;
  case CMR_MATROID_DEC_TYPE_GRAPH:
    fprintf(stream, "graphic matrix with %zu nodes and %zu edges {", CMRgraphNumNodes(dec->graph), CMRgraphNumEdges(dec->graph));
  break;
  case CMR_MATROID_DEC_TYPE_COGRAPH:
    fprintf(stream, "cographic matrix with %zu nodes and %zu edges {", CMRgraphNumNodes(dec->cograph), CMRgraphNumEdges(dec->cograph));
  break;
  case CMR_MATROID_DEC_TYPE_PLANAR:
    assert(CMRgraphNumEdges(dec->graph) == CMRgraphNumEdges(dec->cograph));
    fprintf(stream, "planar matrix with %zu nodes, %zu faces and %zu edges {", CMRgraphNumNodes(dec->graph),
      CMRgraphNumNodes(dec->cograph), CMRgraphNumEdges(dec->graph));
  break;
  case CMR_MATROID_DEC_TYPE_SERIES_PARALLEL:
    if (dec->numChildren)
      fprintf(stream, "matrix with %zu series-parallel reductions; 1 child {", dec->numSeriesParallelReductions);
    else
      fprintf(stream, "series-parallel matrix {");
  break;
  case CMR_MATROID_DEC_TYPE_R10:
    fprintf(stream, "matrix representing R10 {");
  break;
  case CMR_MATROID_DEC_TYPE_FANO:
    fprintf(stream, "matrix representing F_7 {");
  break;
  case CMR_MATROID_DEC_TYPE_FANO_DUAL:
    fprintf(stream, "matrix representing F_7^* {");
  break;
  case CMR_MATROID_DEC_TYPE_K5:
    fprintf(stream, "matrix representing K_5 {");
  break;
  case CMR_MATROID_DEC_TYPE_K5_DUAL:
    fprintf(stream, "matrix representing K_5^* {");
  break;
  case CMR_MATROID_DEC_TYPE_K33:
    fprintf(stream, "matrix representing K_{3,3} {");
  break;
  case CMR_MATROID_DEC_TYPE_K33_DUAL:
    fprintf(stream, "matrix representing K_{3,3}^* {");
  break;
  case CMR_MATROID_DEC_TYPE_SUBMATRIX:
    fprintf(stream, "submatrix {");
  break;
  case CMR_MATROID_DEC_TYPE_PIVOTS:
    fprintf(stream, "pivot {");
  break;
  case CMR_MATROID_DEC_TYPE_DETERMINANT:
    fprintf(stream, "bad determinant {");
  break;
  default:
    fprintf(stream, "[invalid CMR_MATROID_DEC type]");
    fflush(stream);
    return CMR_ERROR_INVALID;
  break;
  }

  bool isFirst = true;
  if (dec->regularity)
  {
    fprintf(stream, (dec->regularity > 0) ? "regular" : "irregular");
    isFirst = false;
  }
  if (dec->graphicness)
  {
    fprintf(stream, "%s%s", isFirst ? "" : ",", (dec->graphicness > 0) ? "graphic" : "not graphic");
    isFirst = false;
  }
  if (dec->cographicness)
  {
    fprintf(stream, "%s%s", isFirst ? "" : ",", (dec->cographicness > 0) ? "cographic" : "not cographic");
    isFirst = false;
  }
  fprintf(stream, "}\n");

  if (printMatrices)
  {
    if (dec->matrix || dec->transpose)
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with matrix:\n\n");
    }
    if (dec->matrix)
      CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stream, '0', false) );
    else if (dec->transpose)
    {
      CMR_CHRMAT* matrix = NULL;
      CMR_CALL( CMRchrmatTranspose(cmr, dec->transpose, &matrix) );
      CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stream, '0', false) );
      CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    }
    else
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "without matrix.\n");
    }
  }

  if (printGraphs)
  {
    if (dec->graph)
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with graph:\n\n");
      CMR_CALL( CMRgraphPrint(stream, dec->graph) );
    }
    if (dec->cograph)
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with cograph:\n\n");
      CMR_CALL( CMRgraphPrint(stream, dec->cograph) );
    }
  }

  if (printReductions && dec->type == CMR_MATROID_DEC_TYPE_SERIES_PARALLEL)
  {
    for (size_t i = 0; i < indent; ++i)
      fputc(' ', stream);
    fprintf(stream, "with series-parallel reductions:\n");
    for (size_t i = 0; i < dec->numSeriesParallelReductions; ++i)
      fprintf(stream, " %s", CMRspReductionString(dec->seriesParallelReductions[i], NULL) );
    fputc('\n', stream);
  }

  for (size_t c = 0; c < dec->numChildren; ++c)
  {
    CMR_CALL( CMRmatroiddecPrint(cmr, dec->children[c], stream, indent + 2, printMatrices, printGraphs,
      printReductions) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRmatroiddecFree(CMR* cmr, CMR_MATROID_DEC** pdec)
{
  assert(cmr);
  assert(pdec);

  CMR_MATROID_DEC* dec = *pdec;
  if (!dec)
    return CMR_OKAY;

  for (size_t c = 0; c < dec->numChildren; ++c)
    CMR_CALL( CMRmatroiddecFree(cmr, &dec->children[c]) );

  CMR_CALL( CMRchrmatFree(cmr, &dec->matrix) );
  CMR_CALL( CMRchrmatFree(cmr, &dec->transpose) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->children) );

  CMR_CALL( CMRfreeBlockArray(cmr, &dec->rowsChild) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->rowsParent) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->rowsRootElement) );

  CMR_CALL( CMRfreeBlockArray(cmr, &dec->columnsChild) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->columnsParent) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->columnsRootElement) );

  CMR_CALL( CMRgraphFree(cmr, &dec->graph) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphForest) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphCoforest) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphArcsReversed) );

  CMR_CALL( CMRgraphFree(cmr, &dec->cograph) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographForest) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographCoforest) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographArcsReversed) );

  CMR_CALL( CMRfreeBlockArray(cmr, &dec->seriesParallelReductions) );

  CMR_CALL( CMRfreeBlockArray(cmr, &dec->pivotRows) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->pivotColumns) );

  CMR_CALL( CMRdensebinmatrixFree(cmr, &dec->denseMatrix) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->denseRowsOriginal) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->denseColumnsOriginal) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsRowsDense) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsColumnsDense) );

  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsSequenceNumRows) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsSequenceNumColumns) );

  CMR_CALL( CMRchrmatFree(cmr, &dec->nestedMinorsMatrix) );
  CMR_CALL( CMRchrmatFree(cmr, &dec->nestedMinorsTranspose) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsRowsOriginal) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsColumnsOriginal) );

  CMR_CALL( CMRfreeBlock(cmr, pdec) );

  return CMR_OKAY;
}

static
CMR_ERROR createNode(
  CMR* cmr,
  CMR_MATROID_DEC** pdec,
  bool isTernary,
  CMR_MATROID_DEC_TYPE type,
  CMR_MATROID_DEC* parent,
  size_t numRows,
  size_t numColumns
)
{
  assert(cmr);
  assert(pdec);

  CMR_CALL( CMRallocBlock(cmr, pdec) );
  CMR_MATROID_DEC* dec = *pdec;
  dec->type = type;
  dec->isTernary = isTernary;
  dec->parent = parent;

  dec->regularity = 0;
  dec->graphicness = 0;
  dec->cographicness = 0;

  dec->matrix = NULL;
  dec->transpose = NULL;

  dec->numChildren = 0;
  dec->children = NULL;

  dec->numRows = numRows;
  dec->rowsChild = NULL;
  dec->rowsParent = NULL;
  dec->rowsRootElement = NULL;
  if (numRows)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &dec->rowsChild, numRows) );
    CMR_CALL( CMRallocBlockArray(cmr, &dec->rowsParent, numRows) );
    CMR_CALL( CMRallocBlockArray(cmr, &dec->rowsRootElement, numRows) );
    for (size_t row = 0; row < numRows; ++row)
    {
      dec->rowsChild[row] = SIZE_MAX;
      dec->rowsParent[row] = SIZE_MAX;
    }
  }

  dec->numColumns = numColumns;
  dec->columnsChild = NULL;
  dec->columnsParent = NULL;
  dec->columnsRootElement = NULL;
  if (numColumns)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &dec->columnsChild, numColumns) );
    CMR_CALL( CMRallocBlockArray(cmr, &dec->columnsParent, numColumns) );
    CMR_CALL( CMRallocBlockArray(cmr, &dec->columnsRootElement, numColumns) );
    for (size_t column = 0; column < numColumns; ++column)
    {
      dec->columnsChild[column] = SIZE_MAX;
      dec->columnsParent[column] = SIZE_MAX;
    }
  }

  dec->testedTwoConnected = false;
  dec->testedR10 = false;

  dec->graph = NULL;
  dec->graphForest = NULL;
  dec->graphCoforest = NULL;
  dec->graphArcsReversed = NULL;

  dec->cograph = NULL;
  dec->cographForest = NULL;
  dec->cographCoforest = NULL;
  dec->cographArcsReversed = NULL;

  dec->testedSeriesParallel = false;
  dec->seriesParallelReductions = NULL;
  dec->numSeriesParallelReductions = 0;

  dec->numPivots = 0;
  dec->pivotRows = NULL;
  dec->pivotColumns = NULL;

  dec->denseMatrix = NULL;
  dec->denseRowsOriginal = NULL;
  dec->denseColumnsOriginal = NULL;
  dec->nestedMinorsRowsDense = NULL;
  dec->nestedMinorsColumnsDense = NULL;

  dec->nestedMinorsLength = 0;
  dec->nestedMinorsSequenceNumRows = NULL;
  dec->nestedMinorsSequenceNumColumns = NULL;

  dec->nestedMinorsMatrix = NULL;
  dec->nestedMinorsTranspose = NULL;
  dec->nestedMinorsRowsOriginal = NULL;
  dec->nestedMinorsColumnsOriginal = NULL;

  dec->nestedMinorsLastGraphic = SIZE_MAX;
  dec->nestedMinorsLastCographic = SIZE_MAX;

  return CMR_OKAY;
}

/**
 * \brief Updates rowsParent and columnsParent of \p dec.
 */

static
CMR_ERROR updateRowsColumnsToParent(
  CMR_MATROID_DEC* dec, /**< Decomposition node. */
  size_t* parentRows,   /**< Array indicating parent row of each row of \p dec. */
  size_t* parentColumns /**< Array indicating parent column of each column of \p dec. */
)
{
  assert(dec);
  assert(parentRows);
  assert(parentColumns);

  for (size_t row = 0; row < dec->numRows; ++row)
  {
    size_t parentRow = parentRows[row];
    dec->rowsParent[row] = parentRow;
  }

  for (size_t column = 0; column < dec->numColumns; ++column)
  {
    size_t parentColumn = parentColumns[column];
    dec->columnsParent[column] = parentColumn;
  }

  return CMR_OKAY;
}

/**
 * \brief Updates rowsChild and columnsChild of the parent.
 */

static
CMR_ERROR updateRowsColumnsFromParent(
  CMR_MATROID_DEC* dec,   /**< Decomposition node. */
  size_t* parentRows,     /**< Array indicating parent row of each row of \p dec. */
  size_t firstRow,        /**< First index of \p parentRows to consider. */
  size_t beyondRow,       /**< Beyond index of \p parentRows to consider. */
  size_t* parentColumns,  /**< Array indicating parent column of each column of \p dec. */
  size_t firstColumn,     /**< First index of \p parentColumns to consider. */
  size_t beyondColumn,    /**< Beyond index of \p parentColumns to consider. */
  size_t childIndex       /**< Child index of \p dec in its parent node. */
)
{
  assert(dec);
  assert(parentRows);
  assert(parentColumns);

  for (size_t row = firstRow; row < beyondRow; ++row)
  {
    size_t parentRow = parentRows[row];
    assert(parentRow != SIZE_MAX);
    dec->parent->rowsChild[parentRow] = childIndex;
  }

  for (size_t column = firstColumn; column < beyondColumn; ++column)
  {
    size_t parentColumn = parentColumns[column];
    assert(parentColumn != SIZE_MAX);
    dec->parent->columnsChild[parentColumn] = childIndex;
  }

  return CMR_OKAY;
}

/**
 * \brief Updates matrix of \p dec from parent using rowsParent and columnsParent arrays.
 */

static
CMR_ERROR updateMatrix(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
)
{
  assert(cmr);
  assert(dec);
  assert(dec->parent);
  assert(dec->parent->matrix);

  CMR_CALL( CMRchrmatFilter(cmr, dec->parent->matrix, dec->numRows, dec->rowsParent, dec->numColumns,
    dec->columnsParent, &dec->matrix) );

  return CMR_OKAY;
}


CMR_ERROR CMRmatroiddecCreateMatrixRoot(CMR* cmr, CMR_MATROID_DEC** pdec, bool isTernary, CMR_CHRMAT* matrix)
{
  assert(cmr);
  assert(pdec);
  assert(matrix);

  CMR_CALL( createNode(cmr, pdec, isTernary, CMR_MATROID_DEC_TYPE_UNKNOWN, NULL, matrix->numRows, matrix->numColumns) );
  CMR_MATROID_DEC* dec = *pdec;

  CMR_CALL( CMRchrmatCopy(cmr, matrix, &dec->matrix) );

  assert(matrix->numRows == dec->numRows);
  assert(matrix->numColumns == dec->numColumns);

  return CMR_OKAY;
}

CMR_ERROR CMRmatroiddecCreateChildFromMatrices(CMR* cmr, CMR_MATROID_DEC* parent, size_t childIndex, CMR_CHRMAT* matrix,
  CMR_CHRMAT* transpose, size_t* rowsParent, size_t* columnsParent)
{
  assert(cmr);
  assert(parent);
  assert(childIndex < parent->numChildren);
  assert(matrix);
  assert(rowsParent);
  assert(columnsParent);

  CMR_CALL( createNode(cmr, &parent->children[childIndex], parent->isTernary, CMR_MATROID_DEC_TYPE_UNKNOWN, parent,
    matrix->numRows, matrix->numColumns) );
  CMR_MATROID_DEC* dec = parent->children[childIndex];
  dec->matrix = matrix;
  dec->transpose = transpose;

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t parentRow = rowsParent[row];
    dec->rowsParent[row] = parentRow;
    if (parentRow != SIZE_MAX)
      parent->rowsChild[parentRow] = childIndex;
  }

  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    size_t parentColumn = columnsParent[column];
    dec->columnsParent[column] = parentColumn;
    if (parentColumn != SIZE_MAX)
      parent->columnsChild[parentColumn] = childIndex;
  }

  return CMR_OKAY;
}


CMR_ERROR CMRmatroiddecUpdateOneSum (CMR* cmr, CMR_MATROID_DEC* dec, size_t numChildren)
{
  assert(cmr);
  assert(dec);
  assert(dec->type == CMR_MATROID_DEC_TYPE_UNKNOWN);
  assert(numChildren >= 2);

  dec->type = CMR_MATROID_DEC_TYPE_ONE_SUM;
  dec->numChildren = numChildren;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->children, numChildren) );
  for (size_t c = 0; c < numChildren; ++c)
    dec->children[c] = NULL;

  return CMR_OKAY;
}


CMR_ERROR CMRmatroiddecInitializeParent(CMR* cmr, CMR_MATROID_DEC* dec, size_t parentsChildIndex,
  size_t* rowsToParentRow, size_t* columnsToParentColumn)
{
  assert(cmr);
  assert(dec);
  assert(dec->numRows > 0);
  assert(dec->numColumns > 0);
  assert(rowsToParentRow);
  assert(columnsToParentColumn);

  for (size_t row = 0; row < dec->numRows; ++row)
  {
    size_t parentRow = rowsToParentRow[row];
    dec->rowsParent[row] = parentRow;
    dec->parent->rowsChild[parentRow] = parentsChildIndex;
    dec->rowsRootElement[row] = dec->parent->rowsRootElement[parentRow];
  }

  for (size_t column = 0; column < dec->numColumns; ++column)
  {
    size_t parentColumn = columnsToParentColumn[column];
    dec->columnsParent[column] = parentColumn;
    dec->parent->columnsChild[parentColumn] = parentsChildIndex;
    dec->columnsRootElement[column] = dec->parent->columnsRootElement[parentColumn];
  }

  return CMR_OKAY;
}

CMR_ERROR CMRmatroiddecUpdateSubmatrix(CMR* cmr, CMR_MATROID_DEC* dec, CMR_SUBMAT* submatrix,
  CMR_MATROID_DEC_TYPE type)
{
  assert(cmr);
  assert(dec);
  assert(submatrix);
  assert(dec->matrix);
  assert(submatrix->numRows <= dec->matrix->numRows);
  assert(submatrix->numColumns <= dec->matrix->numColumns);

  if (submatrix->numRows == dec->matrix->numRows && submatrix->numColumns == dec->matrix->numColumns)
  {
    dec->type = type;
  }
  else
  {
    dec->type = CMR_MATROID_DEC_TYPE_SUBMATRIX;
    dec->numChildren = 1;
    assert(dec->children == NULL);
    CMR_CALL( CMRallocBlockArray(cmr, &dec->children, 1) );
    dec->children[0] = NULL;

    CMR_CHRMAT* childMatrix = NULL;
    CMR_CALL( CMRchrmatZoomSubmat(cmr, dec->matrix, submatrix, &childMatrix) );
    CMR_CALL( createNode(cmr, &dec->children[0], dec->isTernary, type, dec, childMatrix->numRows,
      childMatrix->numColumns) );
    dec->children[0]->matrix = childMatrix;
    CMR_CALL( updateRowsColumnsToParent(dec->children[0], submatrix->rows, submatrix->columns) );
    CMR_CALL( updateRowsColumnsFromParent(dec->children[0], submatrix->rows, 0, childMatrix->numRows,
      submatrix->columns, 0, childMatrix->numColumns, 0) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRmatroiddecUpdateTwoSum(CMR* cmr, CMR_MATROID_DEC* dec, CMR_SEPA* separation)
{
  assert(cmr);
  assert(dec);
  assert(separation);

  size_t numBaseRows[2];
  size_t numBaseColumns[2];
  CMR_CALL( CMRsepaComputeSizes(separation, &numBaseRows[0], &numBaseColumns[0], &numBaseRows[1], &numBaseColumns[1]) );

  dec->type = CMR_MATROID_DEC_TYPE_TWO_SUM;
  dec->numChildren = 2;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->children, 2) );
  for (size_t childIndex = 0; childIndex < 2; ++childIndex)
  {
    size_t numExtraRows = 1 - childIndex;
    size_t numExtraColumns = childIndex;

    dec->children[childIndex] = NULL;
    CMR_CALL( createNode(cmr, &dec->children[childIndex], dec->isTernary, CMR_MATROID_DEC_TYPE_UNKNOWN, dec,
      numBaseRows[childIndex] + numExtraRows, numBaseColumns[childIndex] + numExtraColumns) );
    CMR_MATROID_DEC* child = dec->children[childIndex];

    /* Compute parent rows of child. */
    size_t* parentRows = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &parentRows, numBaseRows[childIndex] + numExtraRows) );
    size_t childRow = 0;

    for (size_t row = 0; row < separation->numRows; ++row)
    {
      if ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == (childIndex == 0 ? CMR_SEPA_FIRST : CMR_SEPA_SECOND))
        parentRows[childRow++] = row;
    }
    assert(childRow == numBaseRows[childIndex]);

    if (childIndex == 0)
    {
      for (size_t row = 0; row < separation->numRows; ++row)
      {
        if (((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          && (separation->rowsFlags[row] & CMR_SEPA_FLAG_RANK1) )
        {
          parentRows[childRow++] = row;
          break;
        }
      }
    }

    /* Compute parent columns of child. */
    size_t* parentColumns = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &parentColumns, numBaseColumns[childIndex] + numExtraColumns) );
    size_t childColumn = 0;

    if (childIndex == 1)
    {
      for (size_t column = 0; column < separation->numColumns; ++column)
      {
        if (((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          && (separation->columnsFlags[column] & CMR_SEPA_FLAG_RANK1) )
        {
          parentColumns[childColumn++] = column;
          break;
        }
      }
    }

    for (size_t column = 0; column < separation->numColumns; ++column)
    {
      if ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == (childIndex == 0 ? CMR_SEPA_FIRST : CMR_SEPA_SECOND))
        parentColumns[childColumn++] = column;
    }
    assert(childColumn == numBaseColumns[childIndex] + numExtraColumns);

    CMR_CALL( updateRowsColumnsToParent(child, parentRows, parentColumns) );
    CMR_CALL( updateRowsColumnsFromParent(child, parentRows, 0, numBaseRows[childIndex], parentColumns, childIndex,
      numBaseColumns[childIndex] + childIndex, childIndex) );

    CMR_CALL( CMRfreeStackArray(cmr, &parentColumns) );
    CMR_CALL( CMRfreeStackArray(cmr, &parentRows) );

    CMR_CALL( updateMatrix(cmr, child) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRmatroiddecUpdatePivots(CMR* cmr, CMR_MATROID_DEC* dec, size_t numPivots, size_t* pivotRows,
  size_t* pivotColumns, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose)
{
  assert(cmr);
  assert(dec);
  assert(numPivots > 0);
  assert(pivotRows);
  assert(pivotColumns);
  assert(matrix);

  dec->type = CMR_MATROID_DEC_TYPE_PIVOTS;
  dec->numChildren = 1;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->children, 1) );
  dec->children[0] = NULL;
  CMR_CALL( createNode(cmr, &dec->children[0], dec->isTernary, CMR_MATROID_DEC_TYPE_UNKNOWN, dec, dec->numRows,
    dec->numColumns) );
  dec->children[0]->matrix = matrix;
  dec->children[0]->transpose = transpose;

  dec->numPivots = numPivots;
  CMR_CALL( CMRduplicateBlockArray(cmr, &dec->pivotRows, numPivots, pivotRows) );
  CMR_CALL( CMRduplicateBlockArray(cmr, &dec->pivotColumns, numPivots, pivotColumns) );

  return CMR_OKAY;
}

CMR_ERROR CMRmatroiddecUpdateThreeSumInit(CMR* cmr, CMR_MATROID_DEC* dec)
{
  assert(cmr);
  assert(dec);

  dec->type = CMR_MATROID_DEC_TYPE_THREE_SUM;
  dec->numChildren = 2;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->children, 2) );
  dec->children[0] = NULL;
  dec->children[1] = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRmatroiddecUpdateThreeSumCreateWideFirstChild(CMR* cmr, CMR_MATROID_DEC* dec, CMR_SEPA* separation,
  size_t* rowsToChild, size_t* columnsToChild, size_t numChildBaseRows, size_t numChildBaseColumns, size_t extraRow,
  size_t extraColumn1, size_t extraColumn2, int8_t extraEntry)
{
  assert(cmr);
  assert(dec);
  assert(dec->matrix);
  assert(dec->transpose);
  assert(separation);
  assert(rowsToChild);
  assert(columnsToChild);
  assert(numChildBaseRows < dec->numRows);
  assert(numChildBaseColumns < dec->numColumns);
  assert(extraRow < dec->numRows);
  assert(extraColumn1 < dec->numColumns);
  assert(extraColumn2 < dec->numColumns);
  assert(extraEntry >= -1 || extraEntry <= 1);

  /* We first create the extra column densely. */
  size_t* parentRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentRows, numChildBaseRows + 1) );
  size_t* parentColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentColumns, numChildBaseColumns + 2) );
  int8_t* denseExtraColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseExtraColumn, numChildBaseRows) );
  for (size_t column = 0; column < numChildBaseRows; ++column)
    denseExtraColumn[column] = 0;
  size_t first = dec->transpose->rowSlice[extraColumn1];
  size_t beyond = dec->transpose->rowSlice[extraColumn1 + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t row = dec->transpose->entryColumns[e];
    size_t childRow = rowsToChild[row];
    if (childRow < numChildBaseRows)
      denseExtraColumn[childRow] = dec->transpose->entryValues[e];
  }

  /* Main matrix. */
  size_t childEntry = 0;
  CMR_CHRMAT* childMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &childMatrix, numChildBaseRows + 1, numChildBaseColumns + 2,
    dec->matrix->numNonzeros + dec->matrix->numRows) );
  for (size_t row = 0; row < dec->numRows; ++row)
  {
    size_t childRow = rowsToChild[row];
    if (childRow >= numChildBaseRows)
      continue;

    parentRows[childRow] = row;
    childMatrix->rowSlice[childRow] = childEntry;

    /* Main entries. */
    size_t first = dec->matrix->rowSlice[row];
    size_t beyond = dec->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = dec->matrix->entryColumns[e];
      size_t childColumn = columnsToChild[column];
      if (childColumn >= numChildBaseColumns)
        continue;

      childMatrix->entryColumns[childEntry] = childColumn;
      childMatrix->entryValues[childEntry] = dec->matrix->entryValues[e];
      ++childEntry;
    }

    /* Dense entries. */
    if (denseExtraColumn[childRow])
    {
      childMatrix->entryColumns[childEntry] = numChildBaseColumns;
      childMatrix->entryValues[childEntry] = denseExtraColumn[childRow];
      ++childEntry;
      childMatrix->entryColumns[childEntry] = numChildBaseColumns + 1;
      childMatrix->entryValues[childEntry] = denseExtraColumn[childRow];
      ++childEntry;
    }
  }
  childMatrix->rowSlice[numChildBaseRows] = childEntry;

  /* Extra row */
  first = dec->matrix->rowSlice[extraRow];
  beyond = dec->matrix->rowSlice[extraRow + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t column = dec->matrix->entryColumns[e];
    size_t childColumn = columnsToChild[column];
    if (childColumn >= numChildBaseColumns)
      continue;

    childMatrix->entryColumns[childEntry] = childColumn;
    childMatrix->entryValues[childEntry] = dec->matrix->entryValues[e];
    ++childEntry;
  }
  parentRows[numChildBaseRows] = extraRow;
  parentColumns[numChildBaseColumns] = extraColumn1;
  parentColumns[numChildBaseColumns + 1] = extraColumn2;
  for (size_t column = 0; column < dec->numColumns; ++column)
  {
    size_t childColumn = columnsToChild[column];
    if (childColumn < numChildBaseColumns)
      parentColumns[childColumn] = column;
  }

  /* Extra entry */
  childMatrix->entryColumns[childEntry] = numChildBaseColumns + 1;
  childMatrix->entryValues[childEntry] = extraEntry;
  ++childEntry;

  /* Finalize */
  childMatrix->rowSlice[numChildBaseRows + 1] = childEntry;
  childMatrix->numNonzeros = childEntry;

  /* Create the actual decomposition node. */
  CMR_CALL( createNode(cmr, &dec->children[0], dec->isTernary, CMR_MATROID_DEC_TYPE_UNKNOWN, dec, childMatrix->numRows,
    childMatrix->numColumns) );
  CMR_MATROID_DEC* child = dec->children[0];
  child->matrix = childMatrix;

  CMR_CALL( updateRowsColumnsToParent(child, parentRows, parentColumns) );
  CMR_CALL( updateRowsColumnsFromParent(child, parentRows, 0, numChildBaseRows, parentColumns, 0, numChildBaseColumns,
    0) );

  CMR_CALL( CMRfreeStackArray(cmr, &denseExtraColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentRows) );

#if defined(CMR_DEBUG)
  CMRdbgMsg(10, "Wide first child matrix:\n");
  CMRchrmatPrintDense(cmr, childMatrix, stdout, '0', true);
  for (size_t childRow = 0; childRow < child->numRows; ++childRow)
    CMRdbgMsg(12, "Child row r%zu corresponds to parent row r%zu.\n", childRow + 1, child->rowsParent[childRow] + 1);
  for (size_t childColumn = 0; childColumn < child->numColumns; ++childColumn)
  {
    CMRdbgMsg(12, "Child column c%zu corresponds to parent column c%zu.\n", childColumn + 1,
      child->columnsParent[childColumn] + 1);
  }
#endif /* CMR_DEBUG */

  return CMR_OKAY;
}

CMR_ERROR CMRmatroiddecUpdateThreeSumCreateWideSecondChild(CMR* cmr, CMR_MATROID_DEC* dec, CMR_SEPA* separation,
  size_t* rowsToChild, size_t* columnsToChild, size_t numChildBaseRows, size_t numChildBaseColumns, size_t extraRow,
  size_t extraColumn1, size_t extraColumn2, int8_t extraEntry)
{
  assert(cmr);
  assert(dec);
  assert(dec->matrix);
  assert(dec->transpose);
  assert(separation);
  assert(rowsToChild);
  assert(columnsToChild);
  assert(numChildBaseRows < dec->numRows);
  assert(numChildBaseColumns < dec->numColumns);
  assert(extraRow < dec->numRows);
  assert(extraColumn1 < dec->numColumns);
  assert(extraColumn2 < dec->numColumns);
  assert(extraEntry >= -1 || extraEntry <= 1);

  /* We first create the extra column densely. */
  size_t* parentRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentRows, numChildBaseRows + 1) );
  size_t* parentColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentColumns, numChildBaseColumns + 2) );
  int8_t* denseExtraColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseExtraColumn, numChildBaseRows) );
  for (size_t column = 0; column < numChildBaseRows; ++column)
    denseExtraColumn[column] = 0;
  size_t first = dec->transpose->rowSlice[extraColumn1];
  size_t beyond = dec->transpose->rowSlice[extraColumn1 + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t row = dec->transpose->entryColumns[e];
    size_t childRow = rowsToChild[row];
    if (childRow < numChildBaseRows)
      denseExtraColumn[childRow] = dec->transpose->entryValues[e];
  }

  /* Main matrix. */
  size_t childEntry = 0;
  CMR_CHRMAT* childMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &childMatrix, numChildBaseRows + 1, numChildBaseColumns + 2,
    dec->matrix->numNonzeros + dec->matrix->numRows) );

  /* Extra entry */
  childMatrix->entryColumns[childEntry] = 0;
  childMatrix->entryValues[childEntry] = extraEntry;
  ++childEntry;

  /* Extra row */
  childMatrix->rowSlice[0] = 0;
  first = dec->matrix->rowSlice[extraRow];
  beyond = dec->matrix->rowSlice[extraRow + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t column = dec->matrix->entryColumns[e];
    size_t childColumn = columnsToChild[column];
    if (childColumn >= numChildBaseColumns)
      continue;

    childMatrix->entryColumns[childEntry] = childColumn + 2;
    childMatrix->entryValues[childEntry] = dec->matrix->entryValues[e];
    ++childEntry;
  }

  for (size_t row = 0; row < dec->numRows; ++row)
  {
    size_t childRow = rowsToChild[row];
    if (childRow >= numChildBaseRows)
      continue;

    parentRows[childRow + 1] = row;
    childMatrix->rowSlice[childRow + 1] = childEntry;

    /* Dense entries. */
    if (denseExtraColumn[childRow])
    {
      childMatrix->entryColumns[childEntry] = 0;
      childMatrix->entryValues[childEntry] = denseExtraColumn[childRow];
      ++childEntry;
      childMatrix->entryColumns[childEntry] = 1;
      childMatrix->entryValues[childEntry] = denseExtraColumn[childRow];
      ++childEntry;
    }

    /* Main entries. */
    size_t first = dec->matrix->rowSlice[row];
    size_t beyond = dec->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = dec->matrix->entryColumns[e];
      size_t childColumn = columnsToChild[column];
      if (childColumn >= numChildBaseColumns)
        continue;

      childMatrix->entryColumns[childEntry] = childColumn + 2;
      childMatrix->entryValues[childEntry] = dec->matrix->entryValues[e];
      ++childEntry;
    }
  }
  childMatrix->rowSlice[numChildBaseRows + 1] = childEntry;

  parentRows[0] = extraRow;
  parentColumns[0] = extraColumn1;
  parentColumns[1] = extraColumn2;
  for (size_t column = 0; column < dec->numColumns; ++column)
  {
    size_t childColumn = columnsToChild[column];
    if (childColumn < numChildBaseColumns)
      parentColumns[childColumn + 2] = column;
  }

  /* Finalize */
  childMatrix->rowSlice[numChildBaseRows + 1] = childEntry;
  childMatrix->numNonzeros = childEntry;

  /* Create the actual decomposition node. */
  CMR_CALL( createNode(cmr, &dec->children[1], dec->isTernary, CMR_MATROID_DEC_TYPE_UNKNOWN, dec, childMatrix->numRows,
    childMatrix->numColumns) );
  CMR_MATROID_DEC* child = dec->children[1];
  child->matrix = childMatrix;

  CMR_CALL( updateRowsColumnsToParent(child, parentRows, parentColumns) );
  CMR_CALL( updateRowsColumnsFromParent(child, parentRows, 1, numChildBaseRows + 1, parentColumns, 2,
    numChildBaseColumns + 2, 1) );

  CMR_CALL( CMRfreeStackArray(cmr, &denseExtraColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentRows) );

#if defined(CMR_DEBUG)
  CMRdbgMsg(10, "Wide second child matrix:\n");
  CMRchrmatPrintDense(cmr, childMatrix, stdout, '0', true);
  for (size_t childRow = 0; childRow < child->numRows; ++childRow)
    CMRdbgMsg(12, "Child row r%zu corresponds to parent row r%zu.\n", childRow + 1, child->rowsParent[childRow] + 1);
  for (size_t childColumn = 0; childColumn < child->numColumns; ++childColumn)
  {
    CMRdbgMsg(12, "Child column c%zu corresponds to parent column c%zu.\n", childColumn + 1,
      child->columnsParent[childColumn] + 1);
  }
#endif /* CMR_DEBUG */

  assert(extraEntry || !"NOT IMPLEMENTED");

  return CMR_OKAY;
}

#define MIN_IF_EXISTS(currentValue, child, childValue) \
  ( (child != NULL) \
    ? ( ((childValue) < (currentValue)) ? (childValue) : (currentValue) ) \
    : ( (0 < (currentValue)) ? 0 : (currentValue) ) \
  )

CMR_ERROR CMRmatroiddecSetAttributes(CMR_MATROID_DEC* dec)
{
  assert(dec);

  CMRdbgMsg(2, "Setting regularity/(co)graphicness attributes for a node of type %d.\n", dec->type);

  /* Recursively set attributes first. */
  for (size_t childIndex = 0; childIndex < dec->numChildren; ++childIndex)
  {
    if (dec->children[childIndex])
    {
      CMR_CALL( CMRmatroiddecSetAttributes(dec->children[childIndex]) );
      CMRdbgMsg(4, "Child %zu of type %d has flags r=%d,g=%d,c=%d.\n", childIndex, dec->children[childIndex]->type,
        dec->children[childIndex]->regularity, dec->children[childIndex]->graphicness,
        dec->children[childIndex]->cographicness);
    }
  }

  switch(dec->type)
  {
  case CMR_MATROID_DEC_TYPE_UNKNOWN:
    dec->regularity = 0;
    dec->graphicness = 0;
    dec->cographicness = 0;
  break;
  case CMR_MATROID_DEC_TYPE_SERIES_PARALLEL:
  case CMR_MATROID_DEC_TYPE_PLANAR:
    dec->regularity = 1;
    dec->graphicness = 1;
    dec->cographicness = 1;
  break;
  case CMR_MATROID_DEC_TYPE_IRREGULAR:
  case CMR_MATROID_DEC_TYPE_FANO:
  case CMR_MATROID_DEC_TYPE_FANO_DUAL:
  case CMR_MATROID_DEC_TYPE_DETERMINANT:
    dec->regularity = -1;
    dec->graphicness = -1;
    dec->cographicness = -1;
  break;
  case CMR_MATROID_DEC_TYPE_PIVOTS:
  case CMR_MATROID_DEC_TYPE_SUBMATRIX:
  case CMR_MATROID_DEC_TYPE_ONE_SUM:
  case CMR_MATROID_DEC_TYPE_TWO_SUM:
  case CMR_MATROID_DEC_TYPE_THREE_SUM:
    if (dec->regularity == 0)
    {
      dec->regularity = 1;
      for (size_t childIndex = 0; childIndex < dec->numChildren; ++childIndex)
      {
        CMR_MATROID_DEC* child = dec->children[childIndex];
        dec->regularity = MIN_IF_EXISTS(dec->regularity, child, child->regularity);
      }
    }
    if (dec->graphicness == 0)
    {
      dec->graphicness = 1;
      for (size_t childIndex = 0; childIndex < dec->numChildren; ++childIndex)
      {
        CMR_MATROID_DEC* child = dec->children[childIndex];
        dec->graphicness = MIN_IF_EXISTS(dec->graphicness, child, child->graphicness);
      }
    }
    if (dec->cographicness == 0)
    {
      dec->cographicness = 1;
      for (size_t childIndex = 0; childIndex < dec->numChildren; ++childIndex)
      {
        CMR_MATROID_DEC* child = dec->children[childIndex];
        dec->cographicness = MIN_IF_EXISTS(dec->cographicness, child, child->cographicness);
      }
    }
  break;
  case CMR_MATROID_DEC_TYPE_GRAPH:
  case CMR_MATROID_DEC_TYPE_K5:
  case CMR_MATROID_DEC_TYPE_K33:
    dec->regularity = 1;
    dec->graphicness = 1;
  break;
  case CMR_MATROID_DEC_TYPE_COGRAPH:
  case CMR_MATROID_DEC_TYPE_K5_DUAL:
  case CMR_MATROID_DEC_TYPE_K33_DUAL:
    dec->regularity = 1;
    dec->cographicness = 1;
  break;
  case CMR_MATROID_DEC_TYPE_R10:
    dec->regularity = 1;
    dec->graphicness = -1;
    dec->cographicness = -1;
  break;
  default:
    assert(0 != "Handling of matroid decomposition type not implemented!");
  }

  return CMR_OKAY;
}

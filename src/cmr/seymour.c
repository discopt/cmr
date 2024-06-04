// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/seymour.h>

#include "env_internal.h"
#include "matroid_internal.h"
#include "matrix_internal.h"
#include "listmatrix.h"

#include <assert.h>
#include <string.h>

bool CMRseymourIsTernary(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->isTernary;
}

bool CMRseymourThreeSumDistributedRanks(CMR_SEYMOUR_NODE* node)
{
  assert(node);
  assert(node->type == CMR_SEYMOUR_NODE_TYPE_THREE_SUM);

  return node->threesumFlags & CMR_SEYMOUR_NODE_THREESUM_FLAG_DISTRIBUTED_RANKS;
}

bool CMRseymourThreeSumConcentratedRank(CMR_SEYMOUR_NODE* node)
{
  assert(node);
  assert(node->type == CMR_SEYMOUR_NODE_TYPE_THREE_SUM);

  return node->threesumFlags & CMR_SEYMOUR_NODE_THREESUM_FLAG_CONCENTRATED_RANK;
}

CMR_CHRMAT* CMRseymourGetMatrix(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->matrix;
}

bool CMRseymourHasTranspose(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->transpose != NULL;
}

CMR_CHRMAT* CMRseymourGetTranspose(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->transpose;
}


size_t CMRseymourNumChildren(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->numChildren;
}

CMR_SEYMOUR_NODE* CMRseymourChild(CMR_SEYMOUR_NODE* node, size_t childIndex)
{
  assert(node);
  assert(childIndex < node->numChildren);

  return node->children[childIndex];
}

CMR_SEYMOUR_NODE_TYPE CMRseymourType(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->type;
}

int8_t CMRseymourGraphicness(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->graphicness;
}

int8_t CMRseymourCographicness(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->cographicness;
}

int8_t CMRseymourRegularity(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->regularity;
}

size_t CMRseymourNumRows(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->numRows;
}

CMR_ELEMENT* CMRseymourChildRowsToParent(CMR_SEYMOUR_NODE* node, size_t childIndex)
{
  assert(node);
  assert(childIndex < node->numChildren);

  return node->childRowsToParent[childIndex];

  return NULL;
}

CMR_ELEMENT* CMRseymourChildColumnsToParent(CMR_SEYMOUR_NODE* node, size_t childIndex)
{
  assert(node);
  assert(childIndex < node->numChildren);

  return node->childColumnsToParent[childIndex];

  return NULL;
}

size_t CMRseymourNumColumns(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->numColumns;
}

CMR_GRAPH* CMRseymourGraph(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->graph;
}

CMR_GRAPH_EDGE* CMRseymourGraphForest(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->graphForest;
}

size_t CMRseymourGraphSizeForest(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  if (node->matrix)
    return node->numRows;
  else if (node->transpose)
    return node->numColumns;
  else if (node->graph)
    return CMRgraphNumNodes(node->graph) - 1;
  else
    return SIZE_MAX;
}

CMR_GRAPH_EDGE* CMRseymourGraphCoforest(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->graphCoforest;
}

size_t CMRseymourGraphSizeCoforest(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  if (node->matrix)
    return node->numColumns;
  else if (node->transpose)
    return node->numRows;
  else if (node->graph)
    return CMRgraphNumEdges(node->graph) + 1 - CMRgraphNumNodes(node->graph);
  else
    return SIZE_MAX;
}

bool* CMRseymourGraphArcsReversed(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->graphArcsReversed;
}

CMR_GRAPH* CMRseymourCograph(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->cograph;
}

size_t CMRseymourCographSizeForest(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  if (node->matrix)
    return node->numColumns;
  else if (node->transpose)
    return node->numRows;
  else if (node->cograph)
    return CMRgraphNumNodes(node->cograph) - 1;
  else
    return SIZE_MAX;
}

CMR_GRAPH_EDGE* CMRseymourCographForest(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->cographForest;
}

size_t CMRseymourCographSizeCoforest(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  if (node->matrix)
    return node->numRows;
  else if (node->transpose)
    return node->numColumns;
  else if (node->cograph)
    return CMRgraphNumEdges(node->cograph) + 1 - CMRgraphNumNodes(node->cograph);
  else
    return SIZE_MAX;
}

CMR_GRAPH_EDGE* CMRseymourCographCoforest(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->cographCoforest;
}

bool* CMRseymourCographArcsReversed(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->cographArcsReversed;
}

size_t CMRseymourNumPivots(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->numPivots;
}

size_t * CMRseymourPivotRows(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->pivotRows;
}

size_t * CMRseymourPivotColumns(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->pivotColumns;
}

CMR_ERROR CMRseymourPrintChild(CMR* cmr, CMR_SEYMOUR_NODE* child, CMR_SEYMOUR_NODE* parent, size_t childIndex,
  FILE* stream, size_t indent, bool printChildren, bool printParentElements, bool printMatrices, bool printGraphs,
  bool printReductions, bool printPivots)
{
  assert(cmr);
  assert(stream);

  /* Indent. */
  for (size_t i = 0; i < indent; ++i)
    fputc(' ', stream);

  if (!child)
  {
    fprintf(stream, "<NULL>\n");
    return CMR_OKAY;
  }

  fprintf(stream, "%zux%zu ", child->numRows, child->numColumns);
  switch (child->type)
  {
  case CMR_SEYMOUR_NODE_TYPE_IRREGULAR:
    fprintf(stream, "irregular node {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_UNKNOWN:
    fprintf(stream, "unknown node {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_ONE_SUM:
    fprintf(stream, "1-sum node with %zu children {", child->numChildren);
  break;
  case CMR_SEYMOUR_NODE_TYPE_TWO_SUM:
    fprintf(stream, "2-sum node {");
    break;
  case CMR_SEYMOUR_NODE_TYPE_THREE_SUM:
    fprintf(stream, "3-sum node {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_GRAPH:
    fprintf(stream, "graphic matrix with %zu nodes and %zu edges {", CMRgraphNumNodes(child->graph),
      CMRgraphNumEdges(child->graph));
  break;
  case CMR_SEYMOUR_NODE_TYPE_COGRAPH:
    fprintf(stream, "cographic matrix with %zu nodes and %zu edges {", CMRgraphNumNodes(child->cograph),
      CMRgraphNumEdges(child->cograph));
  break;
  case CMR_SEYMOUR_NODE_TYPE_PLANAR:
    assert(CMRgraphNumEdges(child->graph) == CMRgraphNumEdges(child->cograph));
    fprintf(stream, "planar matrix with %zu nodes, %zu faces and %zu edges {", CMRgraphNumNodes(child->graph),
      CMRgraphNumNodes(child->cograph), CMRgraphNumEdges(child->graph));
  break;
  case CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL:
    if (child->numChildren)
      fprintf(stream, "matrix with %zu series-parallel reductions; 1 child {", child->numSeriesParallelReductions);
    else
      fprintf(stream, "series-parallel matrix {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_R10:
    fprintf(stream, "matrix representing R10 {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_FANO:
    fprintf(stream, "matrix representing F_7 {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_FANO_DUAL:
    fprintf(stream, "matrix representing F_7^* {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_K5:
    fprintf(stream, "matrix representing K_5 {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_K5_DUAL:
    fprintf(stream, "matrix representing K_5^* {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_K33:
    fprintf(stream, "matrix representing K_{3,3} {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_K33_DUAL:
    fprintf(stream, "matrix representing K_{3,3}^* {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_SUBMATRIX:
    fprintf(stream, "submatrix node {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_PIVOTS:
    fprintf(stream, "pivot node {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_DETERMINANT:
    fprintf(stream, "bad determinant {");
  break;
  default:
    fprintf(stream, "[invalid CMR_MATROID_DEC type]");
    fflush(stream);
    return CMR_ERROR_INVALID;
  break;
  }

  bool isFirst = true;
  if (child->regularity)
  {
    fprintf(stream, (child->regularity > 0) ? "regular" : "irregular");
    isFirst = false;
  }
  if (child->graphicness)
  {
    fprintf(stream, "%s%s", isFirst ? "" : ",", (child->graphicness > 0) ? "graphic" : "not graphic");
    isFirst = false;
  }
  if (child->cographicness)
  {
    fprintf(stream, "%s%s", isFirst ? "" : ",", (child->cographicness > 0) ? "cographic" : "not cographic");
    isFirst = false;
  }
  fprintf(stream, "}\n");

  if (printReductions && child->type == CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL)
  {
    for (size_t i = 0; i < indent; ++i)
      fputc(' ', stream);
    fprintf(stream, "with series-parallel reductions:\n");
    for (size_t i = 0; i < indent; ++i)
      fputc(' ', stream);
    for (size_t i = 0; i < child->numSeriesParallelReductions; ++i)
      fprintf(stream, "%s%s", (i == 0) ? " " : ", ", CMRspReductionString(child->seriesParallelReductions[i], NULL) );
    fputc('\n', stream);
  }

  if (printPivots && child->type == CMR_SEYMOUR_NODE_TYPE_PIVOTS)
  {
    for (size_t i = 0; i < indent; ++i)
      fputc(' ', stream);
    fprintf(stream, "with %zu pivot%s:", child->numPivots, child->numPivots == 1 ? "" : "s");
    for (size_t i = 0; i < child->numPivots; ++i)
      fprintf(stream, "%s r%zu,c%zu", i == 0 ? "" : ",", child->pivotRows[i]+1, child->pivotColumns[i]+1);
    fputc('\n', stream);
  }

  if (printParentElements && parent)
  {
    assert(childIndex < parent->numChildren);

    if (parent->childRowsToParent && parent->childRowsToParent[childIndex])
    {
      CMR_ELEMENT* rowsToParent = parent->childRowsToParent[childIndex];

      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with mapping of rows to parent:");
      for (size_t row = 0; row < child->numRows; ++row)
      {
        if (CMRelementIsRow(rowsToParent[row]))
          fprintf(stream, " r%zu", CMRelementToRowIndex(rowsToParent[row])+1);
        else if (CMRelementIsColumn(rowsToParent[row]))
          fprintf(stream, " c%zu", CMRelementToColumnIndex(rowsToParent[row])+1);
        else
          fprintf(stream, " N/A");
      }
      fprintf(stream, "\n");
    }
    if (parent->childColumnsToParent && parent->childColumnsToParent[childIndex])
    {
      CMR_ELEMENT* columnsToParent = parent->childColumnsToParent[childIndex];

      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with mapping of columns to parent's columns:");
      for (size_t column = 0; column < child->numColumns; ++column)
      {
        if (CMRelementIsRow(columnsToParent[column]))
          fprintf(stream, " r%zu", CMRelementToRowIndex(columnsToParent[column])+1);
        else if (CMRelementIsColumn(columnsToParent[column]))
          fprintf(stream, " c%zu", CMRelementToColumnIndex(columnsToParent[column])+1);
        else
          fprintf(stream, " N/A");
      }
      fprintf(stream, "\n");
    }
  }

  if (printMatrices)
  {
    if (child->matrix || child->transpose)
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with matrix:\n");
    }
    if (child->matrix)
      CMR_CALL( CMRchrmatPrintDense(cmr, child->matrix, stream, '0', false) );
    else if (child->transpose)
    {
      CMR_CHRMAT* matrix = NULL;
      CMR_CALL( CMRchrmatTranspose(cmr, child->transpose, &matrix) );
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
    if (child->graph)
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with graph:\n\n");
      CMR_CALL( CMRgraphPrint(child->graph, stream) );
    }
    if (child->cograph)
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with cograph:\n\n");
      CMR_CALL( CMRgraphPrint(child->cograph, stream) );
    }
  }

  if (printChildren)
  {
    for (size_t c = 0; c < child->numChildren; ++c)
    {
      for (size_t i = 0; i < indent; ++i)
          fputc(' ', stream);
      if (child->numChildren == 1)
        fprintf(stream, "Unique child:\n");
      else
        fprintf(stream, "Child #%zu:\n", c+1);
      CMR_CALL( CMRseymourPrintChild(cmr, child->children[c], child, c, stream, indent + 2, printChildren,
        printParentElements, printMatrices, printGraphs, printReductions, printPivots) );
    }
  }

  return CMR_OKAY;
}

CMR_ERROR CMRseymourPrint(CMR* cmr, CMR_SEYMOUR_NODE* node, FILE* stream, bool printChildren, bool printParentElements,
  bool printMatrices, bool printGraphs, bool printReductions, bool printPivots)
{
  assert(cmr);
  assert(stream);

  CMR_CALL( CMRseymourPrintChild(cmr, node, NULL, SIZE_MAX, stream, 0, printChildren, printParentElements,
    printMatrices, printGraphs, printReductions, printPivots) );

  return CMR_OKAY;
}

CMR_ERROR CMRseymourCapture(CMR* cmr, CMR_SEYMOUR_NODE* node)
{
  assert(cmr);
  assert(node);

  node->used++;

  return CMR_OKAY;
}

CMR_ERROR CMRseymourRelease(CMR* cmr, CMR_SEYMOUR_NODE** pnode)
{
  assert(pnode);

  CMR_SEYMOUR_NODE* node = *pnode;
  assert(node);

  CMRdbgMsg(0, "CMRseymourRelease called for a node with usage %zu.\n", node->used);

  node->used--;
  if (!node->used)
  {
    /* Release recursively. */
    for (size_t c = 0; c < node->numChildren; ++c)
    {
      CMR_CALL( CMRseymourRelease(cmr, &node->children[c]) );
      CMR_CALL( CMRfreeBlockArray(cmr, &node->childRowsToParent[c]) );
      CMR_CALL( CMRfreeBlockArray(cmr, &node->childColumnsToParent[c]) );
    }

    CMR_CALL( CMRchrmatFree(cmr, &node->matrix) );
    CMR_CALL( CMRchrmatFree(cmr, &node->transpose) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->children) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->childRowsToParent) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->childColumnsToParent) );

    CMR_CALL( CMRfreeBlockArray(cmr, &node->rowsToChild) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->columnsToChild) );

    CMR_CALL( CMRgraphFree(cmr, &node->graph) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->graphForest) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->graphCoforest) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->graphArcsReversed) );

    CMR_CALL( CMRgraphFree(cmr, &node->cograph) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->cographForest) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->cographCoforest) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->cographArcsReversed) );

    CMR_CALL( CMRfreeBlockArray(cmr, &node->seriesParallelReductions) );

    CMR_CALL( CMRfreeBlockArray(cmr, &node->pivotRows) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->pivotColumns) );

    CMR_CALL( CMRdensebinmatrixFree(cmr, &node->denseMatrix) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->denseRowsOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->denseColumnsOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->nestedMinorsRowsDense) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->nestedMinorsColumnsDense) );

    CMR_CALL( CMRfreeBlockArray(cmr, &node->nestedMinorsSequenceNumRows) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->nestedMinorsSequenceNumColumns) );

    CMR_CALL( CMRchrmatFree(cmr, &node->nestedMinorsMatrix) );
    CMR_CALL( CMRchrmatFree(cmr, &node->nestedMinorsTranspose) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->nestedMinorsRowsOriginal) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->nestedMinorsColumnsOriginal) );

    CMR_CALL( CMRfreeBlock(cmr, pnode) );
  }
  *pnode = NULL;

  return CMR_OKAY;
}


static
CMR_ERROR createNode(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE** pnode,     /**< Pointer for storing the new Seymour decomposition node. */
  bool isTernary,               /**< Whether the node is ternary. */
  CMR_SEYMOUR_NODE_TYPE type,   /**< Type of the new node. */
  size_t numRows,               /**< Number of rows. */
  size_t numColumns             /**< Number of columns. */
)
{
  assert(cmr);
  assert(pnode);

  CMRdbgMsg(10, "Creating a node for a %zu-by-%zu matrix.\n", numRows, numColumns);

  CMR_CALL( CMRallocBlock(cmr, pnode) );
  CMR_SEYMOUR_NODE* node = *pnode;
  node->type = type;
  node->isTernary = isTernary;
  node->used = 1;

  node->regularity = 0;
  node->graphicness = 0;
  node->cographicness = 0;

  node->matrix = NULL;
  node->transpose = NULL;

  node->numChildren = 0;
  node->children = NULL;
  node->childRowsToParent = NULL;
  node->childColumnsToParent = NULL;

  node->numRows = numRows;
  node->rowsToChild = NULL;
  if (numRows)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &node->rowsToChild, numRows) );
    for (size_t row = 0; row < numRows; ++row)
      node->rowsToChild[row] = SIZE_MAX;
  }

  node->numColumns = numColumns;
  node->columnsToChild = NULL;
  if (numColumns)
  {
    CMR_CALL( CMRallocBlockArray(cmr, &node->columnsToChild, numColumns) );
    for (size_t column = 0; column < numColumns; ++column)
      node->columnsToChild[column] = SIZE_MAX;
  }

  node->testedTwoConnected = false;
  node->testedR10 = false;

  node->graph = NULL;
  node->graphForest = NULL;
  node->graphCoforest = NULL;
  node->graphArcsReversed = NULL;

  node->cograph = NULL;
  node->cographForest = NULL;
  node->cographCoforest = NULL;
  node->cographArcsReversed = NULL;

  node->testedSeriesParallel = false;
  node->seriesParallelReductions = NULL;
  node->numSeriesParallelReductions = 0;

  node->numPivots = 0;
  node->pivotRows = NULL;
  node->pivotColumns = NULL;

  node->denseMatrix = NULL;
  node->denseRowsOriginal = NULL;
  node->denseColumnsOriginal = NULL;
  node->nestedMinorsRowsDense = NULL;
  node->nestedMinorsColumnsDense = NULL;

  node->nestedMinorsLength = 0;
  node->nestedMinorsSequenceNumRows = NULL;
  node->nestedMinorsSequenceNumColumns = NULL;

  node->nestedMinorsMatrix = NULL;
  node->nestedMinorsTranspose = NULL;
  node->nestedMinorsRowsOriginal = NULL;
  node->nestedMinorsColumnsOriginal = NULL;

  node->nestedMinorsLastGraphic = SIZE_MAX;
  node->nestedMinorsLastCographic = SIZE_MAX;

  return CMR_OKAY;
}

CMR_ERROR CMRseymourCloneUnknown(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SEYMOUR_NODE** pclone)
{
  assert(cmr);
  assert(node);
  assert(pclone);

  CMR_CALL( createNode(cmr, pclone, node->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, node->numRows, node->numColumns) );
    CMR_SEYMOUR_NODE* clone = *pclone;

  CMR_CALL( CMRchrmatCopy(cmr, node->matrix, &clone->matrix) );

  return CMR_OKAY;
}

/**
 * \brief Allocates and sets childRowsToParent and childColumnsToParent of the child of \p parent indicated by
 *        \p childIndex.
 */

static
CMR_ERROR updateRowsColumnsToParent(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* parent,   /**< Seymour decomposition parent node. */
  size_t childIndex,          /**< Child index. */
  size_t* parentRows,         /**< Array indicating parent row of each row of the child. */
  size_t* parentColumns       /**< Array indicating parent column of each column of the child. */
)
{
  assert(parent);
  assert(parentRows);
  assert(parentColumns);
  assert(parent->children);
  assert(childIndex < parent->numChildren);

    CMR_SEYMOUR_NODE* child = parent->children[childIndex];
  assert(child);
  assert(parent->childRowsToParent[childIndex] == NULL);
  assert(parent->childColumnsToParent[childIndex] == NULL);

  CMR_CALL( CMRallocBlockArray(cmr, &parent->childRowsToParent[childIndex], child->numRows) );
  CMR_CALL( CMRallocBlockArray(cmr, &parent->childColumnsToParent[childIndex], child->numColumns) );

  for (size_t row = 0; row < child->numRows; ++row)
  {
    size_t parentRow = parentRows[row];
    parent->childRowsToParent[childIndex][row] = parentRow < SIZE_MAX ? CMRrowToElement(parentRow) : 0;
  }

  for (size_t column = 0; column < child->numColumns; ++column)
  {
    size_t parentColumn = parentColumns[column];
    parent->childColumnsToParent[childIndex][column] = parentColumn < SIZE_MAX ? CMRcolumnToElement(parentColumn) : 0;
  }

  return CMR_OKAY;
}

/**
 * \brief Updates rowsToChild and columnsToChild of the parent.
 */

static
CMR_ERROR updateRowsColumnsToChild(
  CMR_SEYMOUR_NODE* parent, /**< Seymour decomposition node of parent. */
  size_t childIndex,        /**< An index of a child of \p parent. */
  size_t* parentRows,       /**< Array indicating parent row of each row of the child. */
  size_t firstRow,          /**< First index of \p parentRows to consider. */
  size_t beyondRow,         /**< Beyond index of \p parentRows to consider. */
  size_t* parentColumns,    /**< Array indicating parent column of each column of the child. */
  size_t firstColumn,       /**< First index of \p parentColumns to consider. */
  size_t beyondColumn       /**< Beyond index of \p parentColumns to consider. */
)
{
  assert(parent);
  assert(parentRows);
  assert(parentColumns);

  for (size_t row = firstRow; row < beyondRow; ++row)
  {
    size_t parentRow = parentRows[row];
    assert(parentRow != SIZE_MAX);
    parent->rowsToChild[parentRow] = childIndex;
  }

  for (size_t column = firstColumn; column < beyondColumn; ++column)
  {
    size_t parentColumn = parentColumns[column];
    assert(parentColumn != SIZE_MAX);
    parent->columnsToChild[parentColumn] = childIndex;
  }

  return CMR_OKAY;
}

/**
 * \brief Updates matrix of \p dec from parent using rowsParent and columnsParent arrays.
 */

static
CMR_ERROR updateChildMatrix(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* parent, /**< Parent node in Seymour decomposition. */
  size_t childIndex         /**< Index of child to update. */
)
{
  assert(cmr);
  assert(parent);
  assert(childIndex < parent->numChildren);

  CMR_SEYMOUR_NODE* child = parent->children[childIndex];

  size_t* rows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rows, child->numRows) );
  for (size_t row = 0; row < child->numRows; ++row)
  {
    assert(CMRelementIsRow(parent->childRowsToParent[childIndex][row]));
    rows[row] = CMRelementToRowIndex(parent->childRowsToParent[childIndex][row]);
  }

  size_t* columns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columns, child->numColumns) );
  for (size_t column = 0; column < child->numColumns; ++column)
  {
    assert(CMRelementIsColumn(parent->childColumnsToParent[childIndex][column]));
    columns[column] = CMRelementToColumnIndex(parent->childColumnsToParent[childIndex][column]);
  }

  CMR_CALL( CMRchrmatFilter(cmr, parent->matrix, child->numRows, rows, child->numColumns, columns, &child->matrix) );

  CMR_CALL( CMRfreeStackArray(cmr, &columns) );
  CMR_CALL( CMRfreeStackArray(cmr, &rows) );

  return CMR_OKAY;
}


CMR_ERROR CMRseymourCreateMatrixRoot(CMR* cmr, CMR_SEYMOUR_NODE** pnode, bool isTernary, CMR_CHRMAT* matrix)
{
  assert(cmr);
  assert(pnode);
  assert(matrix);

  CMR_CALL( createNode(cmr, pnode, isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, matrix->numRows, matrix->numColumns) );
  CMR_SEYMOUR_NODE* node = *pnode;

  CMR_CALL( CMRchrmatCopy(cmr, matrix, &node->matrix) );

  assert(matrix->numRows == node->numRows);
  assert(matrix->numColumns == node->numColumns);

  return CMR_OKAY;
}

CMR_ERROR CMRseymourSetNumChildren(CMR* cmr, CMR_SEYMOUR_NODE* node, size_t numChildren)
{
  assert(cmr);
  assert(node);

  assert(node->numChildren == 0); /* We assume that there were no children so far. */

  node->numChildren = numChildren;
  CMR_CALL( CMRallocBlockArray(cmr, &node->children, numChildren) );
  CMR_CALL( CMRallocBlockArray(cmr, &node->childRowsToParent, numChildren) );
  CMR_CALL( CMRallocBlockArray(cmr, &node->childColumnsToParent, numChildren) );
  for (size_t c = 0; c < numChildren; ++c)
  {
    node->children[c] = NULL;
    node->childRowsToParent[c] = NULL;
    node->childColumnsToParent[c] = NULL;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRseymourCreateChildFromMatrices(CMR* cmr, CMR_SEYMOUR_NODE* parent, size_t childIndex, CMR_CHRMAT* matrix,
  CMR_CHRMAT* transpose, CMR_ELEMENT* rowsToParent, CMR_ELEMENT* columnsToParent)
{
  assert(cmr);
  assert(parent);
  assert(childIndex < parent->numChildren);
  assert(matrix);
  assert(rowsToParent);
  assert(columnsToParent);

  CMR_CALL( createNode(cmr, &parent->children[childIndex], parent->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN,
    matrix->numRows, matrix->numColumns) );
  CMR_SEYMOUR_NODE* node = parent->children[childIndex];
  node->matrix = matrix;
  node->transpose = transpose;

  CMR_CALL( CMRallocBlockArray(cmr, &parent->childRowsToParent[childIndex], node->numRows) );
  CMR_CALL( CMRallocBlockArray(cmr, &parent->childColumnsToParent[childIndex], node->numColumns) );

  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    CMR_ELEMENT parentElement = rowsToParent[row];
    parent->childRowsToParent[childIndex][row] = parentElement;
    if (CMRelementIsRow(parentElement))
      parent->rowsToChild[CMRelementToRowIndex(parentElement)] = childIndex;
    else if (CMRelementIsColumn(parentElement))
      parent->columnsToChild[CMRelementToColumnIndex(parentElement)] = childIndex;
  }

  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    CMR_ELEMENT parentElement = columnsToParent[column];
    parent->childColumnsToParent[childIndex][column] = parentElement;
    if (CMRelementIsRow(parentElement))
      parent->rowsToChild[CMRelementToRowIndex(parentElement)] = childIndex;
    else if (CMRelementIsColumn(parentElement))
      parent->columnsToChild[CMRelementToColumnIndex(parentElement)] = childIndex;
  }

  return CMR_OKAY;
}


CMR_ERROR CMRseymourUpdateOneSum (CMR* cmr, CMR_SEYMOUR_NODE* node, size_t numChildren)
{
  assert(cmr);
  assert(node);
  assert(node->type == CMR_SEYMOUR_NODE_TYPE_UNKNOWN);
  assert(numChildren >= 2);

  node->type = CMR_SEYMOUR_NODE_TYPE_ONE_SUM;

  CMR_CALL( CMRseymourSetNumChildren(cmr, node, numChildren) );

  return CMR_OKAY;
}


// CMR_ERROR CMRseymourInitializeParent(CMR* cmr, CMR_MATROID_DEC* dec, CMR_MATROID_DEC* parent, size_t childIndex,
//   size_t* rowsToParentRow, size_t* columnsToParentColumn)
// {
//   assert(cmr);
//   assert(dec);
//   assert(dec->numRows > 0);
//   assert(dec->numColumns > 0);
//   assert(rowsToParentRow);
//   assert(columnsToParentColumn);
//
//   for (size_t row = 0; row < dec->numRows; ++row)
//   {
//     size_t parentRow = rowsToParentRow[row];
//     dec->rowsParent[row] = CMRrowToElement(parentRow);
//     dec->parent->rowsChild[parentRow] = parentsChildIndex;
//   }
//
//   for (size_t column = 0; column < dec->numColumns; ++column)
//   {
//     size_t parentColumn = columnsToParentColumn[column];
//     dec->columnsParent[column] = CMRcolumnToElement(parentColumn);
//     dec->parent->columnsChild[parentColumn] = parentsChildIndex;
//   }
//
//   return CMR_OKAY;
// }

CMR_ERROR CMRseymourUpdateSubmatrix(CMR* cmr, CMR_SEYMOUR_NODE* dec, CMR_SUBMAT* submatrix,
                                       CMR_SEYMOUR_NODE_TYPE type)
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
    dec->type = CMR_SEYMOUR_NODE_TYPE_SUBMATRIX;
    CMR_CALL( CMRseymourSetNumChildren(cmr, dec, 1) );

    CMR_CHRMAT* childMatrix = NULL;
    CMR_CALL( CMRchrmatZoomSubmat(cmr, dec->matrix, submatrix, &childMatrix) );
    CMR_CALL( createNode(cmr, &dec->children[0], dec->isTernary, type, childMatrix->numRows,
      childMatrix->numColumns) );
    dec->children[0]->matrix = childMatrix;

    CMR_CALL( updateRowsColumnsToParent(cmr, dec, 0, submatrix->rows, submatrix->columns) );
    CMR_CALL( updateRowsColumnsToChild(dec, 0, submatrix->rows, 0, childMatrix->numRows,
      submatrix->columns, 0, childMatrix->numColumns) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdateTwoSum(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SEPA* separation)
{
  assert(cmr);
  assert(node);
  assert(separation);

  size_t numBaseRows[2];
  size_t numBaseColumns[2];
  CMR_CALL( CMRsepaComputeSizes(separation, &numBaseRows[0], &numBaseColumns[0], &numBaseRows[1], &numBaseColumns[1]) );

  node->type = CMR_SEYMOUR_NODE_TYPE_TWO_SUM;
  CMR_CALL( CMRseymourSetNumChildren(cmr, node, 2) );
  for (size_t childIndex = 0; childIndex < 2; ++childIndex)
  {
    size_t numExtraRows = 1 - childIndex;
    size_t numExtraColumns = childIndex;

    CMR_CALL( createNode(cmr, &node->children[childIndex], node->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN,
      numBaseRows[childIndex] + numExtraRows, numBaseColumns[childIndex] + numExtraColumns) );

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

    CMR_CALL( updateRowsColumnsToParent(cmr, node, childIndex, parentRows, parentColumns) );
    CMR_CALL( updateRowsColumnsToChild(node, childIndex, parentRows, 0, numBaseRows[childIndex], parentColumns,
      childIndex, numBaseColumns[childIndex] + childIndex) );

    CMR_CALL( CMRfreeStackArray(cmr, &parentColumns) );
    CMR_CALL( CMRfreeStackArray(cmr, &parentRows) );

    CMR_CALL( updateChildMatrix(cmr, node, childIndex) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdatePivots(CMR* cmr, CMR_SEYMOUR_NODE* node, size_t numPivots, size_t* pivotRows,
  size_t* pivotColumns, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose)
{
  assert(cmr);
  assert(node);
  assert(numPivots > 0);
  assert(pivotRows);
  assert(pivotColumns);
  assert(matrix);

  node->type = CMR_SEYMOUR_NODE_TYPE_PIVOTS;
  CMR_CALL( CMRseymourSetNumChildren(cmr, node, 1) );
  CMR_CALL( createNode(cmr, &node->children[0], node->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, node->numRows,
                         node->numColumns) );
  node->children[0]->matrix = matrix;
  node->children[0]->transpose = transpose;

  node->numPivots = numPivots;
  CMR_CALL( CMRduplicateBlockArray(cmr, &node->pivotRows, numPivots, pivotRows) );
  CMR_CALL( CMRduplicateBlockArray(cmr, &node->pivotColumns, numPivots, pivotColumns) );

  CMR_CALL( CMRallocBlockArray(cmr, &node->childRowsToParent[0], node->numRows) );
  CMR_CALL( CMRallocBlockArray(cmr, &node->childColumnsToParent[0], node->numColumns) );
  for (size_t row = 0; row < node->numRows; ++row)
    node->childRowsToParent[0][row] = CMRrowToElement(row);
  for (size_t column = 0; column < node->numColumns; ++column)
    node->childColumnsToParent[0][column] = CMRcolumnToElement(column);
  for (size_t pivot = 0; pivot < numPivots; ++pivot)
  {
    node->childRowsToParent[0][pivotRows[pivot]] = CMRcolumnToElement(pivotColumns[pivot]);
    node->childColumnsToParent[0][pivotColumns[pivot]] = CMRrowToElement(pivotRows[pivot]);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdateThreeSumInit(CMR* cmr, CMR_SEYMOUR_NODE* node)
{
  assert(cmr);
  assert(node);

  node->type = CMR_SEYMOUR_NODE_TYPE_THREE_SUM;
  CMR_CALL( CMRseymourSetNumChildren(cmr, node, 2) );

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdateThreeSumCreateWideFirstChild(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SEPA* separation,
  size_t* rowsToChild, size_t* columnsToChild, size_t numChildBaseRows, size_t numChildBaseColumns, size_t extraRow,
  size_t extraColumn1, size_t extraColumn2, int8_t extraEntry)
{
  assert(cmr);
  assert(node);
  assert(node->matrix);
  assert(node->transpose);
  assert(separation);
  assert(rowsToChild);
  assert(columnsToChild);
  assert(numChildBaseRows < node->numRows);
  assert(numChildBaseColumns < node->numColumns);
  assert(extraRow < node->numRows);
  assert(extraColumn1 < node->numColumns);
  assert(extraColumn2 < node->numColumns);
  assert(extraEntry == -1 || extraEntry == 1);

  /* We first create the extra column densely. */
  size_t* parentRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentRows, numChildBaseRows + 1) );
  size_t* parentColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentColumns, numChildBaseColumns + 2) );
  int8_t* denseExtraColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseExtraColumn, numChildBaseRows) );
  for (size_t column = 0; column < numChildBaseRows; ++column)
    denseExtraColumn[column] = 0;
  size_t first = node->transpose->rowSlice[extraColumn1];
  size_t beyond = node->transpose->rowSlice[extraColumn1 + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t row = node->transpose->entryColumns[e];
    size_t childRow = rowsToChild[row];
    if (childRow < numChildBaseRows)
      denseExtraColumn[childRow] = node->transpose->entryValues[e];
  }

  /* Main matrix. */
  size_t childEntry = 0;
  CMR_CHRMAT* childMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &childMatrix, numChildBaseRows + 1, numChildBaseColumns + 2,
    node->matrix->numNonzeros + node->matrix->numRows) );
  for (size_t row = 0; row < node->numRows; ++row)
  {
    size_t childRow = rowsToChild[row];
    if (childRow >= numChildBaseRows)
      continue;

    parentRows[childRow] = row;
    childMatrix->rowSlice[childRow] = childEntry;

    /* Main entries. */
    size_t first = node->matrix->rowSlice[row];
    size_t beyond = node->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = node->matrix->entryColumns[e];
      size_t childColumn = columnsToChild[column];
      if (childColumn >= numChildBaseColumns)
        continue;

      childMatrix->entryColumns[childEntry] = childColumn;
      childMatrix->entryValues[childEntry] = node->matrix->entryValues[e];
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
  first = node->matrix->rowSlice[extraRow];
  beyond = node->matrix->rowSlice[extraRow + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t column = node->matrix->entryColumns[e];
    size_t childColumn = columnsToChild[column];
    if (childColumn >= numChildBaseColumns)
      continue;

    childMatrix->entryColumns[childEntry] = childColumn;
    childMatrix->entryValues[childEntry] = node->matrix->entryValues[e];
    ++childEntry;
  }
  parentRows[numChildBaseRows] = extraRow;
  parentColumns[numChildBaseColumns] = extraColumn1;
  parentColumns[numChildBaseColumns + 1] = extraColumn2;
  for (size_t column = 0; column < node->numColumns; ++column)
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
  CMR_CALL( createNode(cmr, &node->children[0], node->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, childMatrix->numRows,
    childMatrix->numColumns) );
  CMR_SEYMOUR_NODE* child = node->children[0];
  child->matrix = childMatrix;

  CMR_CALL( updateRowsColumnsToParent(cmr, node, 0, parentRows, parentColumns) );
  CMR_CALL( updateRowsColumnsToChild(node, 0, parentRows, 0, numChildBaseRows, parentColumns, 0, numChildBaseColumns) );

  CMR_CALL( CMRfreeStackArray(cmr, &denseExtraColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentRows) );

#if defined(CMR_DEBUG)
  CMRdbgMsg(10, "Wide first child matrix:\n");
  CMRchrmatPrintDense(cmr, childMatrix, stdout, '0', true);
  for (size_t childRow = 0; childRow < child->numRows; ++childRow)
  {
    CMRdbgMsg(12, "Child row r%zu corresponds to parent %s.\n", childRow + 1,
      CMRelementString(child->rowsParent[childRow], NULL));
  }
  for (size_t childColumn = 0; childColumn < child->numColumns; ++childColumn)
  {
    CMRdbgMsg(12, "Child column c%zu corresponds to parent %s.\n", childColumn + 1,
      CMRelementString(child->columnsParent[childColumn], NULL));
  }
#endif /* CMR_DEBUG */

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdateThreeSumCreateWideSecondChild(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SEPA* separation,
  size_t* rowsToChild, size_t* columnsToChild, size_t numChildBaseRows, size_t numChildBaseColumns, size_t extraRow,
  size_t extraColumn1, size_t extraColumn2, int8_t extraEntry)
{
  assert(cmr);
  assert(node);
  assert(node->matrix);
  assert(node->transpose);
  assert(separation);
  assert(rowsToChild);
  assert(columnsToChild);
  assert(numChildBaseRows < node->numRows);
  assert(numChildBaseColumns < node->numColumns);
  assert(extraRow < node->numRows);
  assert(extraColumn1 < node->numColumns);
  assert(extraColumn2 < node->numColumns);
  assert(extraEntry == -1 || extraEntry == 1);

  /* We first create the extra column densely. */
  size_t* parentRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentRows, numChildBaseRows + 1) );
  size_t* parentColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentColumns, numChildBaseColumns + 2) );
  int8_t* denseExtraColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseExtraColumn, numChildBaseRows) );
  for (size_t column = 0; column < numChildBaseRows; ++column)
    denseExtraColumn[column] = 0;
  size_t first = node->transpose->rowSlice[extraColumn1];
  size_t beyond = node->transpose->rowSlice[extraColumn1 + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t row = node->transpose->entryColumns[e];
    size_t childRow = rowsToChild[row];
    if (childRow < numChildBaseRows)
      denseExtraColumn[childRow] = node->transpose->entryValues[e];
  }

  /* Main matrix. */
  size_t childEntry = 0;
  CMR_CHRMAT* childMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &childMatrix, numChildBaseRows + 1, numChildBaseColumns + 2,
    node->matrix->numNonzeros + node->matrix->numRows) );

  /* Extra entry */
  childMatrix->entryColumns[childEntry] = 0;
  childMatrix->entryValues[childEntry] = extraEntry;
  ++childEntry;

  /* Extra row */
  childMatrix->rowSlice[0] = 0;
  first = node->matrix->rowSlice[extraRow];
  beyond = node->matrix->rowSlice[extraRow + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t column = node->matrix->entryColumns[e];
    size_t childColumn = columnsToChild[column];
    if (childColumn >= numChildBaseColumns)
      continue;

    childMatrix->entryColumns[childEntry] = childColumn + 2;
    childMatrix->entryValues[childEntry] = node->matrix->entryValues[e];
    ++childEntry;
  }

  for (size_t row = 0; row < node->numRows; ++row)
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
    size_t first = node->matrix->rowSlice[row];
    size_t beyond = node->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = node->matrix->entryColumns[e];
      size_t childColumn = columnsToChild[column];
      if (childColumn >= numChildBaseColumns)
        continue;

      childMatrix->entryColumns[childEntry] = childColumn + 2;
      childMatrix->entryValues[childEntry] = node->matrix->entryValues[e];
      ++childEntry;
    }
  }
  childMatrix->rowSlice[numChildBaseRows + 1] = childEntry;

  parentRows[0] = extraRow;
  parentColumns[0] = extraColumn1;
  parentColumns[1] = extraColumn2;
  for (size_t column = 0; column < node->numColumns; ++column)
  {
    size_t childColumn = columnsToChild[column];
    if (childColumn < numChildBaseColumns)
      parentColumns[childColumn + 2] = column;
  }

  /* Finalize */
  childMatrix->rowSlice[numChildBaseRows + 1] = childEntry;
  childMatrix->numNonzeros = childEntry;

  /* Create the actual decomposition node. */
  CMR_CALL( createNode(cmr, &node->children[1], node->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, childMatrix->numRows,
    childMatrix->numColumns) );
  CMR_SEYMOUR_NODE* child = node->children[1];
  child->matrix = childMatrix;

  CMR_CALL( updateRowsColumnsToParent(cmr, node, 1, parentRows, parentColumns) );
  CMR_CALL( updateRowsColumnsToChild(node, 1, parentRows, 1, numChildBaseRows + 1, parentColumns, 2,
    numChildBaseColumns + 2) );

  CMR_CALL( CMRfreeStackArray(cmr, &denseExtraColumn) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentRows) );

#if defined(CMR_DEBUG)
  CMRdbgMsg(10, "Wide second child matrix:\n");
  CMRchrmatPrintDense(cmr, childMatrix, stdout, '0', true);
  for (size_t childRow = 0; childRow < child->numRows; ++childRow)
  {
    CMRdbgMsg(12, "Child row r%zu corresponds to parent %s.\n", childRow + 1,
      CMRelementString(child->rowsParent[childRow], NULL));
  }
  for (size_t childColumn = 0; childColumn < child->numColumns; ++childColumn)
  {
    CMRdbgMsg(12, "Child column c%zu corresponds to parent %s.\n", childColumn + 1,
      CMRelementString(child->columnsParent[childColumn], NULL));
  }
#endif /* CMR_DEBUG */

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdateThreeSumCreateMixedFirstChild(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SEPA* separation,
  size_t* rowsToChild, size_t* columnsToChild, size_t numChildBaseRows, size_t numChildBaseColumns, size_t extraRow1,
  size_t extraRow2, int8_t extraEntry)
{
  assert(cmr);
  assert(node);
  assert(node->matrix);
  assert(node->transpose);
  assert(separation);
  assert(rowsToChild);
  assert(columnsToChild);
  assert(numChildBaseRows < node->numRows);
  assert(numChildBaseColumns < node->numColumns);
  assert(extraRow1 < node->numRows);
  assert(extraRow2 < node->numRows);
  assert(extraEntry == -1 || extraEntry == 1);

  /* Mapping from child rows/columns to parent rows/columns. */
  size_t* parentRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentRows, numChildBaseRows + 2) );
  size_t* parentColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentColumns, numChildBaseColumns + 1) );

  /* Main matrix. */
  size_t childEntry = 0;
  CMR_CHRMAT* childMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &childMatrix, numChildBaseRows + 2, numChildBaseColumns + 1,
    node->matrix->numNonzeros) );
  for (size_t row = 0; row < node->numRows; ++row)
  {
    size_t childRow = rowsToChild[row];
    if (childRow >= numChildBaseRows)
      continue;

    parentRows[childRow] = row;
    childMatrix->rowSlice[childRow] = childEntry;

    /* Main entries. */
    size_t first = node->matrix->rowSlice[row];
    size_t beyond = node->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = node->matrix->entryColumns[e];
      size_t childColumn = columnsToChild[column];
      if (childColumn >= numChildBaseColumns)
        continue;

      childMatrix->entryColumns[childEntry] = childColumn;
      childMatrix->entryValues[childEntry] = node->matrix->entryValues[e];
      ++childEntry;
    }
  }
  childMatrix->rowSlice[numChildBaseRows] = childEntry;

  /* Extra row 1 with a +1 at the end. */
  size_t first = node->matrix->rowSlice[extraRow1];
  size_t beyond = node->matrix->rowSlice[extraRow1 + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t column = node->matrix->entryColumns[e];
    size_t childColumn = columnsToChild[column];
    if (childColumn >= numChildBaseColumns)
      continue;

    childMatrix->entryColumns[childEntry] = childColumn;
    childMatrix->entryValues[childEntry] = node->matrix->entryValues[e];
    ++childEntry;
  }
  childMatrix->entryColumns[childEntry] = numChildBaseColumns;
  childMatrix->entryValues[childEntry] = 1;
  ++childEntry;
  parentRows[numChildBaseRows] = extraRow1;
  childMatrix->rowSlice[numChildBaseRows + 1] = childEntry;

  /* Extra row 2 with extra entry at the end. */
  first = node->matrix->rowSlice[extraRow2];
  beyond = node->matrix->rowSlice[extraRow2 + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t column = node->matrix->entryColumns[e];
    size_t childColumn = columnsToChild[column];
    if (childColumn >= numChildBaseColumns)
      continue;

    childMatrix->entryColumns[childEntry] = childColumn;
    childMatrix->entryValues[childEntry] = node->matrix->entryValues[e];
    ++childEntry;
  }
  childMatrix->entryColumns[childEntry] = numChildBaseColumns;
  childMatrix->entryValues[childEntry] =  extraEntry;
  ++childEntry;
  parentRows[numChildBaseRows + 1] = extraRow2;
  for (size_t column = 0; column < node->numColumns; ++column)
  {
    size_t childColumn = columnsToChild[column];
    if (childColumn < numChildBaseColumns)
      parentColumns[childColumn] = column;
  }
  parentColumns[numChildBaseColumns] = SIZE_MAX;

  /* Finalize */
  childMatrix->rowSlice[numChildBaseRows + 2] = childEntry;
  childMatrix->numNonzeros = childEntry;

  /* Create the actual decomposition node. */
  CMR_CALL( createNode(cmr, &node->children[0], node->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, childMatrix->numRows,
    childMatrix->numColumns) );
  CMR_SEYMOUR_NODE* child = node->children[0];
  child->matrix = childMatrix;

  CMR_CALL( updateRowsColumnsToParent(cmr, node, 0, parentRows, parentColumns) );
  CMR_CALL( updateRowsColumnsToChild(node, 0, parentRows, 0, numChildBaseRows, parentColumns, 0, numChildBaseColumns) );

  CMR_CALL( CMRfreeStackArray(cmr, &parentColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentRows) );

#if defined(CMR_DEBUG)
  CMRdbgMsg(10, "Mixed first child matrix:\n");
  CMRchrmatPrintDense(cmr, childMatrix, stdout, '0', true);
  for (size_t childRow = 0; childRow < child->numRows; ++childRow)
  {
    CMRdbgMsg(12, "Child row r%zu corresponds to parent %s.\n", childRow + 1,
      CMRelementString(child->rowsParent[childRow], NULL));
  }
  for (size_t childColumn = 0; childColumn < child->numColumns; ++childColumn)
  {
    CMRdbgMsg(12, "Child column c%zu corresponds to parent %s.\n", childColumn + 1,
      CMRelementString(child->columnsParent[childColumn], NULL));
  }
#endif /* CMR_DEBUG */

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdateThreeSumCreateMixedSecondChild(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SEPA* separation,
  size_t* rowsToChild, size_t* columnsToChild, size_t numChildBaseRows, size_t numChildBaseColumns, size_t extraColumn1,
  size_t extraColumn2, int8_t extraEntry)
{
  assert(cmr);
  assert(node);
  assert(node->matrix);
  assert(node->transpose);
  assert(separation);
  assert(rowsToChild);
  assert(columnsToChild);
  assert(numChildBaseRows < node->numRows);
  assert(numChildBaseColumns < node->numColumns);
  assert(extraColumn1 < node->numColumns);
  assert(extraColumn2 < node->numColumns);
  assert(extraEntry == -1 || extraEntry == 1);

  /* We first create the two extra columns densely. */
  size_t* parentRows = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentRows, numChildBaseRows + 1) );
  size_t* parentColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &parentColumns, numChildBaseColumns + 2) );

  int8_t* denseExtraColumn1 = NULL;
  int8_t* denseExtraColumn2 = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseExtraColumn1, numChildBaseRows) );
  CMR_CALL( CMRallocStackArray(cmr, &denseExtraColumn2, numChildBaseRows) );
  for (size_t column = 0; column < numChildBaseRows; ++column)
  {
    denseExtraColumn1[column] = 0;
    denseExtraColumn2[column] = 0;
  }
  size_t first = node->transpose->rowSlice[extraColumn1];
  size_t beyond = node->transpose->rowSlice[extraColumn1 + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t row = node->transpose->entryColumns[e];
    size_t childRow = rowsToChild[row];
    if (childRow < numChildBaseRows)
      denseExtraColumn1[childRow] = node->transpose->entryValues[e];
  }
  first = node->transpose->rowSlice[extraColumn2];
  beyond = node->transpose->rowSlice[extraColumn2 + 1];
  for (size_t e = first; e < beyond; ++e)
  {
    size_t row = node->transpose->entryColumns[e];
    size_t childRow = rowsToChild[row];
    if (childRow < numChildBaseRows)
      denseExtraColumn2[childRow] = node->transpose->entryValues[e];
  }

  /* Main matrix. */
  size_t childEntry = 0;
  CMR_CHRMAT* childMatrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &childMatrix, numChildBaseRows + 1, numChildBaseColumns + 2,
    node->matrix->numNonzeros) );

  /* Only entries in first row are the extra entry and a +1. */
  childMatrix->rowSlice[0] = 0;
  childMatrix->entryColumns[childEntry] = 0;
  childMatrix->entryValues[childEntry] = extraEntry;
  ++childEntry;
  childMatrix->entryColumns[childEntry] = 1;
  childMatrix->entryValues[childEntry] = +1;
  ++childEntry;

  for (size_t row = 0; row < node->numRows; ++row)
  {
    size_t childRow = rowsToChild[row];
    if (childRow >= numChildBaseRows)
      continue;

    parentRows[childRow + 1] = row;
    childMatrix->rowSlice[childRow + 1] = childEntry;

    /* Dense entries. */
    if (denseExtraColumn1[childRow])
    {
      childMatrix->entryColumns[childEntry] = 0;
      childMatrix->entryValues[childEntry] = denseExtraColumn1[childRow];
      ++childEntry;
    }
    if (denseExtraColumn2[childRow])
    {
      childMatrix->entryColumns[childEntry] = 1;
      childMatrix->entryValues[childEntry] = denseExtraColumn2[childRow];
      ++childEntry;
    }

    /* Main entries. */
    size_t first = node->matrix->rowSlice[row];
    size_t beyond = node->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = node->matrix->entryColumns[e];
      size_t childColumn = columnsToChild[column];
      if (childColumn >= numChildBaseColumns)
        continue;

      childMatrix->entryColumns[childEntry] = childColumn + 2;
      childMatrix->entryValues[childEntry] = node->matrix->entryValues[e];
      ++childEntry;
    }
  }
  childMatrix->rowSlice[numChildBaseRows + 1] = childEntry;

  parentRows[0] = SIZE_MAX;
  parentColumns[0] = extraColumn1;
  parentColumns[1] = extraColumn2;
  for (size_t column = 0; column < node->numColumns; ++column)
  {
    size_t childColumn = columnsToChild[column];
    if (childColumn < numChildBaseColumns)
      parentColumns[childColumn + 2] = column;
  }

  /* Finalize */
  childMatrix->rowSlice[numChildBaseRows + 1] = childEntry;
  childMatrix->numNonzeros = childEntry;

  /* Create the actual decomposition node. */
  CMR_CALL( createNode(cmr, &node->children[1], node->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, childMatrix->numRows,
    childMatrix->numColumns) );
  CMR_SEYMOUR_NODE* child = node->children[1];
  child->matrix = childMatrix;

  CMR_CALL( updateRowsColumnsToParent(cmr, node, 1, parentRows, parentColumns) );
  CMR_CALL( updateRowsColumnsToChild(node, 1, parentRows, 1, numChildBaseRows + 1, parentColumns, 2,
    numChildBaseColumns + 2) );

  CMR_CALL( CMRfreeStackArray(cmr, &denseExtraColumn2) );
  CMR_CALL( CMRfreeStackArray(cmr, &denseExtraColumn1) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &parentRows) );

#if defined(CMR_DEBUG)
  CMRdbgMsg(10, "Mixed second child matrix:\n");
  CMRchrmatPrintDense(cmr, childMatrix, stdout, '0', true);
  for (size_t childRow = 0; childRow < child->numRows; ++childRow)
  {
    CMRdbgMsg(12, "Child row r%zu corresponds to parent %s.\n", childRow + 1,
      CMRelementString(child->rowsParent[childRow], NULL));
  }
  for (size_t childColumn = 0; childColumn < child->numColumns; ++childColumn)
  {
    CMRdbgMsg(12, "Child column c%zu corresponds to parent %s.\n", childColumn + 1,
      CMRelementString(child->columnsParent[childColumn], NULL));
  }
#endif /* CMR_DEBUG */

  return CMR_OKAY;
}

#define MIN_IF_EXISTS(currentValue, child, childValue) \
  ( (child != NULL) \
    ? ( ((childValue) < (currentValue)) ? (childValue) : (currentValue) ) \
    : ( (0 < (currentValue)) ? 0 : (currentValue) ) \
  )

CMR_ERROR CMRseymourSetAttributes(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  CMRdbgMsg(2, "Setting regularity/(co)graphicness attributes for a node of type %d.\n", node->type);

  /* Recursively set attributes first. */
  for (size_t childIndex = 0; childIndex < node->numChildren; ++childIndex)
  {
    if (node->children[childIndex])
    {
      CMR_CALL( CMRseymourSetAttributes(node->children[childIndex]) );
      CMRdbgMsg(4, "Child %zu of type %d has flags r=%d,g=%d,c=%d.\n", childIndex, node->children[childIndex]->type,
        node->children[childIndex]->regularity, node->children[childIndex]->graphicness,
        node->children[childIndex]->cographicness);
    }
  }

  switch(node->type)
  {
  case CMR_SEYMOUR_NODE_TYPE_UNKNOWN:
    node->regularity = 0;
    /* We do not set (co)graphicness because its absense may have been determined without setting the type. */
  break;
  case CMR_SEYMOUR_NODE_TYPE_PLANAR:
    node->regularity = 1;
    node->graphicness = 1;
    node->cographicness = 1;
  break;
  case CMR_SEYMOUR_NODE_TYPE_IRREGULAR:
  case CMR_SEYMOUR_NODE_TYPE_FANO:
  case CMR_SEYMOUR_NODE_TYPE_FANO_DUAL:
  case CMR_SEYMOUR_NODE_TYPE_DETERMINANT:
    node->regularity = -1;
    node->graphicness = -1;
    node->cographicness = -1;
  break;
  case CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL:
    if (node->numChildren)
    {
      node->regularity = node->children[0]->regularity;
      node->graphicness = node->children[0]->graphicness;
      node->cographicness = node->children[0]->cographicness;
    }
    else
    {
      node->regularity = 1;
      node->graphicness = 1;
      node->cographicness = 1;
    }
    break;
  case CMR_SEYMOUR_NODE_TYPE_PIVOTS:
  case CMR_SEYMOUR_NODE_TYPE_SUBMATRIX:
  case CMR_SEYMOUR_NODE_TYPE_ONE_SUM:
  case CMR_SEYMOUR_NODE_TYPE_TWO_SUM:
  case CMR_SEYMOUR_NODE_TYPE_THREE_SUM:
    if (node->regularity == 0)
    {
      node->regularity = 1;
      for (size_t childIndex = 0; childIndex < node->numChildren; ++childIndex)
      {
        CMR_SEYMOUR_NODE* child = node->children[childIndex];
        node->regularity = MIN_IF_EXISTS(node->regularity, child, child->regularity);
      }
    }
    if (node->graphicness == 0)
    {
      node->graphicness = 1;
      for (size_t childIndex = 0; childIndex < node->numChildren; ++childIndex)
      {
        CMR_SEYMOUR_NODE* child = node->children[childIndex];
        node->graphicness = MIN_IF_EXISTS(node->graphicness, child, child->graphicness);
      }
    }
    if (node->cographicness == 0)
    {
      node->cographicness = 1;
      for (size_t childIndex = 0; childIndex < node->numChildren; ++childIndex)
      {
        CMR_SEYMOUR_NODE* child = node->children[childIndex];
        node->cographicness = MIN_IF_EXISTS(node->cographicness, child, child->cographicness);
      }
    }
  break;
  case CMR_SEYMOUR_NODE_TYPE_GRAPH:
  case CMR_SEYMOUR_NODE_TYPE_K5:
  case CMR_SEYMOUR_NODE_TYPE_K33:
    node->regularity = 1;
    node->graphicness = 1;
  break;
  case CMR_SEYMOUR_NODE_TYPE_COGRAPH:
  case CMR_SEYMOUR_NODE_TYPE_K5_DUAL:
  case CMR_SEYMOUR_NODE_TYPE_K33_DUAL:
    node->regularity = 1;
    node->cographicness = 1;
  break;
  case CMR_SEYMOUR_NODE_TYPE_R10:
    node->regularity = 1;
    node->graphicness = -1;
    node->cographicness = -1;
  break;
  default:
    assert(0 != "Handling of matroid decomposition type not implemented!");
  }

  CMRdbgMsg(4, "Finalizing attributes to r=%d,g=%d,c=%d.\n", node->regularity, node->graphicness, node->cographicness);

  return CMR_OKAY;
}


typedef struct
{
    CMR_SEYMOUR_NODE* origin;
    CMR_SEYMOUR_NODE* clone;
} ClonePair;

static
CMR_ERROR cloneRecursively(
  CMR* cmr,                                   /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node,                     /**< Seymour decomposition node. */
  CMR_SEYMOUR_NODE** pclone,                  /**< Pointer for storing the cloned Seymour decomposition node. */
  CMR_LISTHASHTABLE* nodesToClonesHashtable,  /**< Hashtable for the mapping from nodes to cloned nodes. */
  ClonePair** pclonePairs,                    /**< Pointer to array of cloned pairs. */
  size_t* pmemClonePairs,                     /**< Allocated size of \p *pclonePairs. */
  size_t* pnumClonePairs                      /**< Length of \p *pclonePairs. */
)
{
  assert(cmr);
  assert(node);
  assert(pclone);
  assert(nodesToClonesHashtable);
  assert(pclonePairs);
  assert(pmemClonePairs);
  assert(pnumClonePairs);

  assert(*pnumClonePairs <= *pmemClonePairs);

  ClonePair* clonePairs = *pclonePairs;
  assert(clonePairs);

  /* Check if node was already cloned. */
  CMR_LISTHASHTABLE_HASH hash = (CMR_LISTHASHTABLE_HASH) node;
  for (CMR_LISTHASHTABLE_ENTRY entry = CMRlisthashtableFindFirst(nodesToClonesHashtable, hash);
    entry != SIZE_MAX; entry = CMRlisthashtableFindNext(nodesToClonesHashtable, hash, entry))
  {
    size_t clonePairIndex = CMRlisthashtableValue(nodesToClonesHashtable, entry);
    if (clonePairs[clonePairIndex].origin == node)
    {
      *pclone = clonePairs[clonePairIndex].clone;
      return CMR_OKAY;
    }
  }

  /* It was not, so we create a clone. */
    CMR_SEYMOUR_NODE* clone = NULL;
  CMR_CALL( createNode(cmr, &clone, node->isTernary, node->type, node->numRows, node->numColumns) );
  CMR_CALL( CMRchrmatCopy(cmr, node->matrix, &clone->matrix) );
  if (node->graph)
    CMR_CALL( CMRgraphCopy(cmr, node->graph, &clone->graph) );
  if (node->graphForest)
    CMR_CALL( CMRduplicateBlockArray(cmr, &clone->graphForest, CMRgraphNumNodes(node->graph) - 1, node->graphForest) );
  if (node->graphCoforest)
  {
    CMR_CALL( CMRduplicateBlockArray(cmr, &clone->graphCoforest,
      CMRgraphNumEdges(node->graph) - CMRgraphNumNodes(node->graph) + 1, node->graphCoforest) );
  }
  if (node->graphArcsReversed)
  {
    CMR_CALL( CMRduplicateBlockArray(cmr, &clone->graphArcsReversed, CMRgraphMemEdges(node->graph),
                                         node->graphArcsReversed) );
  }
  if (node->cograph)
    CMR_CALL( CMRgraphCopy(cmr, node->graph, &clone->graph) );
  if (node->cographForest)
  {
    CMR_CALL( CMRduplicateBlockArray(cmr, &clone->cographForest, CMRgraphNumNodes(node->cograph) - 1,
                                         node->cographForest) );
  }
  if (node->cographCoforest)
  {
    CMR_CALL( CMRduplicateBlockArray(cmr, &clone->cographCoforest,
      CMRgraphNumEdges(node->cograph) - CMRgraphNumNodes(node->cograph) + 1, node->cographCoforest) );
  }
  if (node->cographArcsReversed)
  {
    CMR_CALL( CMRduplicateBlockArray(cmr, &clone->cographArcsReversed, CMRgraphMemEdges(node->cograph),
                                         node->cographArcsReversed) );
  }

  if (pclone)
    *pclone = clone;

  /* Possibly enlarge the list of clone pairs. */
  if (*pnumClonePairs == *pmemClonePairs)
  {
    *pmemClonePairs *= 2;
    CMR_CALL( CMRreallocBlockArray(cmr, pclonePairs, *pmemClonePairs) );
    clonePairs = *pclonePairs;
  }

  /* Add it to the hash table and to the list of clone pairs. */
  CMR_CALL( CMRlisthashtableInsert(cmr, nodesToClonesHashtable, hash, *pnumClonePairs, NULL) );
  assert(*pnumClonePairs < *pmemClonePairs);
  clonePairs[*pnumClonePairs].origin = node;
  clonePairs[*pnumClonePairs].clone = clone;
  ++(*pnumClonePairs);

  /* Recursively treat the children. */
  CMR_CALL( CMRseymourSetNumChildren(cmr, clone, node->numChildren) );
  for (size_t c = 0; c < node->numChildren; ++c)
  {
    /* Clone the child. */
    CMR_CALL( cloneRecursively(cmr, node->children[c], NULL, nodesToClonesHashtable, pclonePairs, pmemClonePairs,
      pnumClonePairs) );

    /* Then add linking data. */
  }

  return CMR_OKAY;
}

CMR_ERROR CMRregularityCloneSubtrees(CMR* cmr, size_t numSubtrees, CMR_SEYMOUR_NODE** subtreeRoots,
  CMR_SEYMOUR_NODE** clonedSubtrees)
{
  assert(cmr);
  assert(subtreeRoots);

  CMR_LISTHASHTABLE* nodesToClonesHashtable = NULL;
  CMR_CALL( CMRlisthashtableCreate(cmr, &nodesToClonesHashtable, 1024, 256) );

  ClonePair* clonePairs = NULL;
  size_t memClonePairs = (numSubtrees > 32) ? 2 * numSubtrees : 32;
  CMR_CALL( CMRallocBlockArray(cmr, &clonePairs, memClonePairs) );
  size_t numClonePairs = 0;

  /* We call the recursive cloning function. */
  for (size_t i = 0; i < numSubtrees; ++i)
  {
    CMR_CALL( cloneRecursively(cmr, subtreeRoots[i], &clonedSubtrees[i], nodesToClonesHashtable, &clonePairs,
      &memClonePairs, &numClonePairs) );
  }

  CMR_CALL( CMRfreeBlockArray(cmr, &clonePairs) );
  CMR_CALL( CMRlisthashtableFree(cmr, &nodesToClonesHashtable) );

  return CMR_OKAY;
}

// #define CMR_DEBUG /* Uncomment to debug this file. */
// #define CMR_DEBUG_MATRICES /* Uncomment to print actual matrices. */

#include <cmr/seymour.h>

#include "env_internal.h"
#include "seymour_internal.h"
#include "matrix_internal.h"
#include "listmatrix.h"

#include <assert.h>
#include <string.h>
#include <time.h>

CMR_ERROR CMRseymourParamsInit(CMR_SEYMOUR_PARAMS* params)
{
  assert(params);

  params->stopWhenIrregular = false;
  params->stopWhenNongraphic = false;
  params->stopWhenNoncographic = false;
  params->stopWhenNeitherGraphicNorCoGraphic = false;
  params->seriesParallel = true;
  params->planarityCheck = false;
  params->directGraphicness = true;
  params->preferGraphicness = true;
  params->decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_DELTASUM
    | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM;
  params->constructLeafGraphs = false;
  params->constructAllGraphs = false;

  return CMR_OKAY;
}

CMR_ERROR CMRseymourStatsInit(CMR_SEYMOUR_STATS* stats)
{
  assert(stats);

  stats->totalCount = 0;
  stats->totalTime = 0.0;
  CMR_CALL( CMRspStatsInit(&stats->seriesParallel) );
  CMR_CALL( CMRgraphicStatsInit(&stats->graphic) );
  CMR_CALL( CMRnetworkStatsInit(&stats->network) );
  stats->sequenceExtensionCount = 0;
  stats->sequenceExtensionTime = 0.0;
  stats->sequenceGraphicCount = 0;
  stats->sequenceGraphicTime = 0.0;
  stats->enumerationCount = 0;
  stats->enumerationTime = 0.0;
  stats->enumerationCandidatesCount = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRseymourStatsPrint(FILE* stream, CMR_SEYMOUR_STATS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  if (!prefix)
  {
    fprintf(stream, "Seymour decomposition:\n");
    prefix = "  ";
  }

  char subPrefix[256];
  snprintf(subPrefix, 256, "%sseries-parallel ", prefix);
  CMR_CALL( CMRspStatsPrint(stream, &stats->seriesParallel, subPrefix) );
  snprintf(subPrefix, 256, "%s(co)graphic ", prefix);
  CMR_CALL( CMRgraphicStatsPrint(stream, &stats->graphic, subPrefix) );
  snprintf(subPrefix, 256, "%s(co)network ", prefix);
  CMR_CALL( CMRnetworkStatsPrint(stream, &stats->network, subPrefix) );

  fprintf(stream, "%ssequence extensions: %lu in %f seconds\n", prefix, (unsigned long)stats->sequenceExtensionCount,
    stats->sequenceExtensionTime);
  fprintf(stream, "%ssequence (co)graphic: %lu in %f seconds\n", prefix, (unsigned long)stats->sequenceGraphicCount,
    stats->sequenceGraphicTime);
  fprintf(stream, "%senumeration: %lu in %f seconds\n", prefix, (unsigned long)stats->enumerationCount,
    stats->enumerationTime);
  fprintf(stream, "%s3-separation candidates: %lu (%.1fk per second)\n", prefix,
    (unsigned long)stats->enumerationCandidatesCount,
    stats->enumerationTime > 0.0 ? (stats->enumerationCandidatesCount / 1000.0 / stats->enumerationTime) : 0.0);
  fprintf(stream, "%stotal: %lu in %f seconds\n", prefix, (unsigned long)stats->totalCount, stats->totalTime);

  return CMR_OKAY;
}

bool CMRseymourIsTernary(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->isTernary;
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
}

CMR_ELEMENT* CMRseymourChildColumnsToParent(CMR_SEYMOUR_NODE* node, size_t childIndex)
{
  assert(node);
  assert(childIndex < node->numChildren);

  return node->childColumnsToParent[childIndex];
}

size_t* CMRseymourChildSpecialRows(CMR_SEYMOUR_NODE* node, size_t childIndex)
{
  assert(node);
  assert(childIndex < node->numChildren);

  return node->childSpecialRows[childIndex];
}

size_t* CMRseymourChildSpecialColumns(CMR_SEYMOUR_NODE* node, size_t childIndex)
{
  assert(node);
  assert(childIndex < node->numChildren);

  return node->childSpecialColumns[childIndex];
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

size_t CMRseymourGetUsed(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->used;
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
  case CMR_SEYMOUR_NODE_TYPE_ONESUM:
    fprintf(stream, "1-sum node with %zu children {", child->numChildren);
  break;
  case CMR_SEYMOUR_NODE_TYPE_TWOSUM:
    fprintf(stream, "2-sum node {");
    break;
  case CMR_SEYMOUR_NODE_TYPE_DELTASUM:
    fprintf(stream, "delta-sum node {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_THREESUM:
    fprintf(stream, "3-sum node {");
  break;
  case CMR_SEYMOUR_NODE_TYPE_YSUM:
    fprintf(stream, "Y-sum node {");
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
  case CMR_SEYMOUR_NODE_TYPE_PIVOTS:
    fprintf(stream, "pivot node {");
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
    if (parent->childSpecialRows && parent->childSpecialRows[childIndex])
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      if (parent->type == CMR_SEYMOUR_NODE_TYPE_DELTASUM)
      {
        fprintf(stream, "with special rows: r%zu\n", parent->childSpecialRows[childIndex][0]+1);
      }
      else if (parent->type == CMR_SEYMOUR_NODE_TYPE_YSUM)
      {
        fprintf(stream, "with special rows: r%zu, r%zu\n", parent->childSpecialRows[childIndex][0]+1,
          parent->childSpecialRows[childIndex][1]+1);
      }
      else if (parent->type == CMR_SEYMOUR_NODE_TYPE_THREESUM)
      {
        fprintf(stream, "with special rows: r%zu, r%zu, r%zu\n", parent->childSpecialRows[childIndex][0]+1,
          parent->childSpecialRows[childIndex][1]+1, parent->childSpecialRows[childIndex][2]+1);
      }
    }
    if (parent->childColumnsToParent && parent->childColumnsToParent[childIndex])
    {
      CMR_ELEMENT* columnsToParent = parent->childColumnsToParent[childIndex];

      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "with mapping of columns to parent:");
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
    if (parent->childSpecialRows && parent->childSpecialRows[childIndex])
    {
      for (size_t i = 0; i < indent; ++i)
        fputc(' ', stream);
      if (parent->type == CMR_SEYMOUR_NODE_TYPE_DELTASUM)
      {
        fprintf(stream, "with special columns: c%zu, c%zu\n", parent->childSpecialColumns[childIndex][0]+1,
          parent->childSpecialColumns[childIndex][1]+1);
      }
      else if (parent->type == CMR_SEYMOUR_NODE_TYPE_YSUM)
      {
        fprintf(stream, "with special columns: c%zu\n", parent->childSpecialColumns[childIndex][0]+1);
      }
      else if (parent->type == CMR_SEYMOUR_NODE_TYPE_THREESUM)
      {
        fprintf(stream, "with special columns: c%zu, c%zu, c%zu\n", parent->childSpecialColumns[childIndex][0]+1,
          parent->childSpecialColumns[childIndex][1]+1, parent->childSpecialColumns[childIndex][2]+1);
      }
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
  CMR_UNUSED(cmr);

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
      CMR_CALL( CMRfreeBlockArray(cmr, &node->childSpecialRows[c]) );
      CMR_CALL( CMRfreeBlockArray(cmr, &node->childSpecialColumns[c]) );
    }

    CMR_CALL( CMRchrmatFree(cmr, &node->matrix) );
    CMR_CALL( CMRchrmatFree(cmr, &node->transpose) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->children) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->childRowsToParent) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->childColumnsToParent) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->childSpecialRows) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->childSpecialColumns) );

    CMR_CALL( CMRfreeBlockArray(cmr, &node->rowsToChild) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->columnsToChild) );

    for (size_t m = 0; m < node->numMinors; ++m)
      CMR_CALL( CMRminorFree(cmr, &node->minors[m]) );
    CMR_CALL( CMRfreeBlockArray(cmr, &node->minors) );

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
  node->childSpecialRows = NULL;
  node->childSpecialColumns = NULL;

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

  node->memMinors = 0;
  node->numMinors = 0;
  node->minors = NULL;

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

  /* Create submatrix is on the stack. */
  CMR_SUBMAT submat;
  submat.numRows = child->numRows;
  submat.numColumns = child->numColumns;
  submat.rows = rows;
  submat.columns = columns;
  CMR_CALL( CMRchrmatSlice(cmr, parent->matrix, &submat, &child->matrix) );

  CMR_CALL( CMRfreeStackArray(cmr, &columns) );
  CMR_CALL( CMRfreeStackArray(cmr, &rows) );

  return CMR_OKAY;
}


CMR_ERROR CMRseymourCreate(CMR* cmr, CMR_SEYMOUR_NODE** pnode, bool isTernary, size_t numRows, size_t numColumns)
{
  assert(cmr);
  assert(pnode);

  CMR_CALL( createNode(cmr, pnode, isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, numRows, numColumns) );

  return CMR_OKAY;
}

CMR_ERROR CMRseymourAddMinor(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_MINOR* minor)
{
  assert(cmr);
  assert(node);
  assert(minor);

  if (node->numMinors == node->memMinors)
  {
    node->memMinors = node->memMinors ? 2 * node->memMinors : 4;
    CMR_CALL( CMRreallocBlockArray(cmr, &node->minors, node->memMinors) );
  }

  node->minors[node->numMinors++] = minor;

  return CMR_OKAY;
}

size_t CMRseymourNumMinors(CMR_SEYMOUR_NODE* node)
{
  assert(node);

  return node->numMinors;
}

CMR_MINOR* CMRseymourMinor(CMR_SEYMOUR_NODE* node, size_t minorIndex)
{
  assert(node);
  assert(minorIndex < node->numMinors);

  return node->minors[minorIndex];
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
  CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialRows, numChildren) );
  CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialColumns, numChildren) );
  for (size_t c = 0; c < numChildren; ++c)
  {
    node->children[c] = NULL;
    node->childRowsToParent[c] = NULL;
    node->childColumnsToParent[c] = NULL;
    node->childSpecialRows[c] = NULL;
    node->childSpecialColumns[c] = NULL;
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


CMR_ERROR CMRseymourUpdateOnesum (CMR* cmr, CMR_SEYMOUR_NODE* node, size_t numChildren)
{
  assert(cmr);
  assert(node);
  assert(node->type == CMR_SEYMOUR_NODE_TYPE_UNKNOWN);
  assert(numChildren >= 2);

  node->type = CMR_SEYMOUR_NODE_TYPE_ONESUM;

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

CMR_ERROR CMRseymourUpdateViolator(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SUBMAT* violator)
{
  assert(cmr);
  assert(node);
  assert(violator);
  assert(node->matrix);
  assert(violator->numRows <= node->matrix->numRows);
  assert(violator->numColumns <= node->matrix->numColumns);

  node->type = CMR_SEYMOUR_NODE_TYPE_IRREGULAR;
  CMR_MINOR* minor = NULL;
  CMR_CALL( CMRminorCreate(cmr, &minor, 0, violator, CMR_MINOR_TYPE_DETERMINANT) );
  CMR_CALL( CMRseymourAddMinor(cmr, node, minor) );

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdateSeriesParallel(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SUBMAT* reducedSubmatrix)
{
  assert(cmr);
  assert(node);
  assert(node->matrix);

  node->type = CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL;
  CMR_CALL( CMRseymourSetNumChildren(cmr, node, 1) );

  CMR_CHRMAT* childMatrix = NULL;
  CMR_CALL( CMRchrmatSlice(cmr, node->matrix, reducedSubmatrix, &childMatrix) );
  CMR_CALL( createNode(cmr, &node->children[0], node->isTernary, CMR_SEYMOUR_NODE_TYPE_UNKNOWN, childMatrix->numRows,
    childMatrix->numColumns) );
  node->children[0]->matrix = childMatrix;

  CMR_CALL( updateRowsColumnsToParent(cmr, node, 0, reducedSubmatrix->rows, reducedSubmatrix->columns) );
  CMR_CALL( updateRowsColumnsToChild(node, 0, reducedSubmatrix->rows, 0, childMatrix->numRows,
    reducedSubmatrix->columns, 0, childMatrix->numColumns) );

  return CMR_OKAY;
}

CMR_ERROR CMRseymourUpdateTwosum(CMR* cmr, CMR_SEYMOUR_NODE* node, CMR_SEPA* separation)
{
  assert(cmr);
  assert(node);
  assert(separation);

  size_t numBaseRows[2];
  size_t numBaseColumns[2];
  CMR_CALL( CMRsepaComputeSizes(separation, &numBaseRows[0], &numBaseColumns[0], &numBaseRows[1], &numBaseColumns[1]) );

  node->type = CMR_SEYMOUR_NODE_TYPE_TWOSUM;
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
  case CMR_SEYMOUR_NODE_TYPE_ONESUM:
  case CMR_SEYMOUR_NODE_TYPE_TWOSUM:
  case CMR_SEYMOUR_NODE_TYPE_THREESUM:
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
  case CMR_SEYMOUR_NODE_TYPE_DELTASUM:
  case CMR_SEYMOUR_NODE_TYPE_YSUM:
    if (node->regularity == 0)
    {
      node->regularity = 1;
      for (size_t childIndex = 0; childIndex < node->numChildren; ++childIndex)
      {
        CMR_SEYMOUR_NODE* child = node->children[childIndex];
        node->regularity = MIN_IF_EXISTS(node->regularity, child, child->regularity);
      }
    }
  break;
  case CMR_SEYMOUR_NODE_TYPE_GRAPH:
    node->regularity = 1;
    node->graphicness = 1;
  break;
  case CMR_SEYMOUR_NODE_TYPE_COGRAPH:
    node->regularity = 1;
    node->cographicness = 1;
  break;
  case CMR_SEYMOUR_NODE_TYPE_R10:
    node->regularity = 1;
    node->graphicness = -1;
    node->cographicness = -1;
  break;
  default:
    assert(0 != "Handling of Seymour decomposition type not implemented!");
  }

  CMRdbgMsg(4, "Finalizing attributes to r=%d,g=%d,c=%d.\n", node->regularity, node->graphicness, node->cographicness);

  return CMR_OKAY;
}

/**
 * \brief Pair of a given Seymour decomposition node and its clone.
 */

typedef struct
{
  CMR_SEYMOUR_NODE* origin;
  CMR_SEYMOUR_NODE* clone;
} ClonePair;

/**
 * \brief Clones \p node to \p *pclone, updates cloning data structures and calls itself recursively.
 *
 * Inserts it into \p nodesToClonesHashtable, and updates the array \p *pclonePairs with the cloned pairs.
 */

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

CMR_ERROR CMRseymourCloneSubtrees(CMR* cmr, size_t numSubtrees, CMR_SEYMOUR_NODE** subtreeRoots,
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


CMR_ERROR CMRregularityTaskCreateRoot(CMR* cmr, CMR_SEYMOUR_NODE* node, DecompositionTask** ptask,
  CMR_SEYMOUR_PARAMS* params, CMR_SEYMOUR_STATS* stats, clock_t startClock, double timeLimit)
{
  assert(cmr);
  assert(node);
  assert(ptask);
  assert(params);

  CMR_CALL( CMRallocBlock(cmr, ptask) );
  DecompositionTask* task = *ptask;

  task->node = node;
  task->next = NULL;

  task->params = params;
  task->stats = stats;
  task->startClock = startClock;
  task->timeLimit = timeLimit;

  return CMR_OKAY;
}

CMR_ERROR CMRregularityTaskFree(CMR* cmr, DecompositionTask** ptask)
{
  assert(cmr);
  assert(ptask);

  if (*ptask == NULL)
    return CMR_OKAY;

  CMR_CALL( CMRfreeBlock(cmr, ptask) );

  return CMR_OKAY;
}

CMR_ERROR CMRregularityQueueCreate(CMR* cmr, DecompositionQueue** pqueue)
{
  assert(cmr);
  assert(pqueue);

  CMR_CALL( CMRallocBlock(cmr, pqueue) );
  DecompositionQueue* queue = *pqueue;
  queue->head = NULL;
  queue->foundIrregularity = false;
  queue->foundNongraphicness = false;
  queue->foundNoncographicness = false;

  return CMR_OKAY;
}

CMR_ERROR CMRregularityQueueFree(CMR* cmr, DecompositionQueue** pqueue)
{
  assert(cmr);
  assert(pqueue);

  DecompositionQueue* queue = *pqueue;
  if (queue == NULL)
    return CMR_OKAY;

  while (queue->head)
  {
    DecompositionTask* task = queue->head;
    queue->head = task->next;
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }

  CMR_CALL( CMRfreeBlock(cmr, pqueue) );

  return CMR_OKAY;
}

bool CMRregularityQueueEmpty(DecompositionQueue* queue)
{
  assert(queue);

  return queue->head == NULL;
}

DecompositionTask* CMRregularityQueueRemove(DecompositionQueue* queue)
{
  assert(queue);

  DecompositionTask* task = queue->head;
  queue->head = task->next;
  task->next = NULL;
  return task;
}

void CMRregularityQueueAdd(DecompositionQueue* queue, DecompositionTask* task)
{
  assert(queue);

  task->next = queue->head;
  queue->head = task;
}

/**
 * \brief Runs a task for processing the associated decomposition node.
 */

static
CMR_ERROR CMRregularityTaskRun(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed tasks. */
)
{
  assert(cmr);
  assert(task);
  assert(queue);

  CMR_ERROR error;

  CMRdbgMsg(2, "Processing task %p.\n", task);

  if (!task->node->testedTwoConnected)
  {
    CMRdbgMsg(4, "Searching for 1-separations.\n");
    error = CMRregularitySearchOnesum(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }
  else if (!task->node->graphicness
    && (task->params->directGraphicness || task->node->matrix->numRows <= 3 || task->node->matrix->numColumns <= 3))
  {
    CMRdbgMsg(4, "Testing directly for %s.\n", task->node->isTernary ? "being network" : "graphicness");
    error = CMRregularityTestGraphicness(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }
  else if (!task->node->cographicness
    && (task->params->directGraphicness || task->node->matrix->numRows <= 3 || task->node->matrix->numColumns <= 3))
  {
    CMRdbgMsg(4, "Testing directly for %s.\n", task->node->isTernary ? "being conetwork" : "cographicness");
    error = CMRregularityTestCographicness(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }
  else if (!task->node->testedR10)
  {
    CMRdbgMsg(4, "Testing for being R_10.\n");
    error = CMRregularityTestR10(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }
  else if (!task->node->testedSeriesParallel)
  {
    CMRdbgMsg(4, "Testing for series-parallel reductions.\n");
    error = CMRregularityDecomposeSeriesParallel(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }
  else if (task->node->denseMatrix)
  {
    CMRdbgMsg(4, "Attempting to construct a sequence of nested minors.\n");
    error = CMRregularityExtendNestedMinorSequence(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }
  else if (task->node->nestedMinorsMatrix && (task->node->nestedMinorsLastGraphic == SIZE_MAX))
  {
    CMRdbgMsg(4, "Testing along the sequence for %s.\n", task->node->isTernary ? "being network" : "graphicness");
    error = CMRregularityNestedMinorSequenceGraphicness(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }
  else if (task->node->nestedMinorsMatrix && (task->node->nestedMinorsLastCographic == SIZE_MAX))
  {
    CMRdbgMsg(4, "Testing along the sequence for %s.\n", task->node->isTernary ? "being conetwork" : "cographicness");
    error = CMRregularityNestedMinorSequenceCographicness(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }
  else
  {
    CMRdbgMsg(4, "Searching for 3-separations along the sequence.\n");
    error = CMRregularityNestedMinorSequenceSearchThreeSeparation(cmr, task, queue);
    if (error != CMR_OKAY && error != CMR_ERROR_TIMEOUT)
      CMR_CALL( error );
  }

  return error;
}

CMR_ERROR CMRseymourDecompose(CMR* cmr, CMR_CHRMAT* matrix, bool ternary, CMR_SEYMOUR_NODE** proot,
  CMR_SEYMOUR_PARAMS* params, CMR_SEYMOUR_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(proot);
  assert(matrix);
  assert(params);

#if defined(CMR_DEBUG_MATRICES)
  CMRdbgMsg(0, "Testing a %s %zux%zu matrix for regularity.\n", ternary ? "ternary" : "binary", matrix->numRows,
    matrix->numColumns);
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
#endif /* CMR_DEBUG_MATRICES */

  CMR_CALL( CMRseymourCreate(cmr, proot, ternary, matrix->numRows, matrix->numColumns) );
  assert(*proot);
  CMR_CALL( CMRchrmatCopy(cmr, matrix, &(*proot)->matrix) );

  CMR_ERROR error = CMRregularityCompleteDecomposition(cmr, *proot, params, stats, timeLimit);
  if (error == CMR_ERROR_TIMEOUT)
  {
    CMR_CALL( CMRseymourRelease(cmr, proot) );
    *proot = NULL;
    return error;
  }
  CMR_CALL( error );

  return CMR_OKAY;
}

CMR_ERROR CMRregularityCompleteDecomposition(CMR* cmr, CMR_SEYMOUR_NODE* subtree, CMR_SEYMOUR_PARAMS* params,
  CMR_SEYMOUR_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(subtree);
  assert(params);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "Completing decomposition tree for a %s %zux%zu matrix.\n",
    subtree->isTernary ? "ternary" : "binary", subtree->matrix->numRows, subtree->matrix->numColumns);
  CMRdbgMsg(0, "Considered subtree belongs to the %zux%zu matrix.\n",
    subtree->matrix->numRows, subtree->matrix->numColumns);
#endif /* CMR_DEBUG */
#if defined(CMR_DEBUG_MATRICES)
  CMR_CALL( CMRchrmatPrintDense(cmr, subtree->matrix, stdout, '0', false) );
#endif /* CMR_DEBUG_MATRICES */

  if (params->decomposeStrategy ==
    (CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_PIVOT | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT))
  {
    return CMR_ERROR_PARAMS;
  }

  clock_t time = clock();
  if (stats)
    stats->totalCount++;

  for (size_t c = 0; c < subtree->numChildren; ++c)
  {
    CMR_CALL( CMRseymourRelease(cmr, &subtree->children[c]) );
    CMR_CALL( CMRfreeBlockArray(cmr, &subtree->childRowsToParent[c]) );
    CMR_CALL( CMRfreeBlockArray(cmr, &subtree->childColumnsToParent[c]) );
    CMR_CALL( CMRfreeBlockArray(cmr, &subtree->childSpecialRows[c]) );
    CMR_CALL( CMRfreeBlockArray(cmr, &subtree->childSpecialColumns[c]) );
  }

  subtree->type = CMR_SEYMOUR_NODE_TYPE_UNKNOWN;

  DecompositionQueue* queue = NULL;
  CMR_CALL( CMRregularityQueueCreate(cmr, &queue) );
  DecompositionTask* decTask = NULL;
  CMR_CALL( CMRregularityTaskCreateRoot(cmr, subtree, &decTask, params, stats, time, timeLimit) );
  CMRregularityQueueAdd(queue, decTask);

  CMRdbgMsg(2, "Main decomposition task is %p\n", decTask);

  CMR_ERROR error = CMR_OKAY;
  while (!CMRregularityQueueEmpty(queue))
  {
    DecompositionTask* task = CMRregularityQueueRemove(queue);

    CMRdbgMsg(2, "Dequeuing task %p\n", task);

    if (params->stopWhenIrregular && queue->foundIrregularity)
    {
      CMRdbgMsg(2, "Clearing task queue due to an irregular node.\n");
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      continue;
    }
    else if (params->stopWhenNongraphic && queue->foundNongraphicness)
    {
      CMRdbgMsg(2, "Clearing task queue due to a nongraphic node.\n");
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      continue;
    }
    else if (params->stopWhenNoncographic && queue->foundNoncographicness)
    {
      CMRdbgMsg(2, "Clearing task queue due to a noncographic node.\n");
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      continue;
    }
    else if (params->stopWhenNeitherGraphicNorCoGraphic && queue->foundNongraphicness && queue->foundNoncographicness)
    {
      CMRdbgMsg(2, "Clearing task queue due to a nongraphic node and a noncographic node.\n");
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      continue;
    }

    error = CMRregularityTaskRun(cmr, task, queue);
    if (error == CMR_ERROR_TIMEOUT)
    {
      CMRdbgMsg(2, "Timeout -> removing task %p.\n", task);
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      break;
    }
    CMR_CALL( error );
  }

  CMR_CALL( CMRregularityQueueFree(cmr, &queue) );

  CMR_CALL( CMRseymourSetAttributes(subtree) );

  if (stats)
    stats->totalTime += (clock() - time) * 1.0 / CLOCKS_PER_SEC;

  return error;
}

CMR_ERROR CMRregularityRefineDecomposition(CMR* cmr, size_t numNodes, CMR_SEYMOUR_NODE** nodes,
  CMR_SEYMOUR_PARAMS* params, CMR_SEYMOUR_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(nodes);
  assert(params);

  CMRdbgMsg(0, "Refining decomposition trees of %zu %s matrices.\n", numNodes,
    nodes[0]->isTernary ? "ternary" : "binary");

  clock_t time = clock();
  if (stats)
    stats->totalCount++;

  DecompositionQueue* queue = NULL;
  CMR_CALL( CMRregularityQueueCreate(cmr, &queue) );
  DecompositionTask* decTask = NULL;

  for (size_t i = 0; i < numNodes; ++i)
  {
    CMR_SEYMOUR_NODE* subtree = nodes[i];
    for (size_t c = 0; c < subtree->numChildren; ++c)
    {
      CMR_CALL( CMRseymourRelease(cmr, &subtree->children[c]) );
      CMR_CALL( CMRfreeBlockArray(cmr, &subtree->childRowsToParent[c]) );
      CMR_CALL( CMRfreeBlockArray(cmr, &subtree->childColumnsToParent[c]) );
      CMR_CALL( CMRfreeBlockArray(cmr, &subtree->childSpecialRows[c]) );
      CMR_CALL( CMRfreeBlockArray(cmr, &subtree->childSpecialColumns[c]) );
    }

    subtree->type = CMR_SEYMOUR_NODE_TYPE_UNKNOWN;

    CMR_CALL( CMRregularityTaskCreateRoot(cmr, subtree, &decTask, params, stats, time, timeLimit) );
    CMRregularityQueueAdd(queue, decTask);
  }

  while (!CMRregularityQueueEmpty(queue))
  {
    DecompositionTask* task = CMRregularityQueueRemove(queue);

    if (params->stopWhenIrregular && queue->foundIrregularity)
    {
      CMRdbgMsg(2, "Clearing task queue due to an irregular node.\n");
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      continue;
    }
    else if (params->stopWhenNongraphic && queue->foundNongraphicness)
    {
      CMRdbgMsg(2, "Clearing task queue due to a nongraphic node.\n");
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      continue;
    }
    else if (params->stopWhenNoncographic && queue->foundNoncographicness)
    {
      CMRdbgMsg(2, "Clearing task queue due to a noncographic node.\n");
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      continue;
    }
    else if (params->stopWhenNeitherGraphicNorCoGraphic && queue->foundNongraphicness && queue->foundNoncographicness)
    {
      CMRdbgMsg(2, "Clearing task queue due to a nongraphic node and a noncographic node.\n");
      CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      continue;
    }

    CMR_CALL( CMRregularityTaskRun(cmr, task, queue) );
  }

  CMR_CALL( CMRregularityQueueFree(cmr, &queue) );

  for (size_t i = 0; i < numNodes; ++i)
    CMR_CALL( CMRseymourSetAttributes(nodes[i]) );

  if (stats)
    stats->totalTime += (clock() - time) * 1.0 / CLOCKS_PER_SEC;

  return CMR_OKAY;
}


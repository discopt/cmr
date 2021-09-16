// #define CMR_DEBUG /* Uncomment to debug this file. */

#include "env_internal.h"
#include "dec_internal.h"
#include "matrix_internal.h"

bool CMRdecHasMatrix(CMR_DEC* dec)
{
  assert(dec);

  return dec->matrix;
}

bool CMRdecHasTranspose(CMR_DEC* dec)
{
  assert(dec);

  return dec->transpose;
}

CMR_CHRMAT* CMRdecGetMatrix(CMR_DEC* dec)
{
  assert(dec);

  return dec->matrix;
}

CMR_CHRMAT* CMRdecGetTranspose(CMR_DEC* dec)
{
  assert(dec);

  return dec->transpose;
}

size_t CMRdecNumChildren(CMR_DEC* dec)
{
  assert(dec);

  return dec->numChildren;
}

CMR_DEC* CMRdecChild(CMR_DEC* dec, size_t childIndex)
{
  assert(dec);
  assert(childIndex < dec->numChildren);

  return dec->children[childIndex];
}

CMR_ERROR CMRdecFree(CMR* cmr, CMR_DEC** pdec)
{
  assert(cmr);
  assert(pdec);

  CMR_DEC* dec = *pdec;
  if (!dec)
    return CMR_OKAY;

  for (size_t c = 0; c < dec->numChildren; ++c)
    CMR_CALL( CMRdecFree(cmr, &dec->children[c]) );

  CMR_CALL( CMRfreeBlockArray(cmr, &dec->children) );
  CMR_CALL( CMRchrmatFree(cmr, &dec->matrix) );
  CMR_CALL( CMRchrmatFree(cmr, &dec->transpose) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->rowsParent) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->columnsParent) );
  CMR_CALL( CMRgraphFree(cmr, &dec->graph) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphForest) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphCoforest) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphArcsReversed) );
  CMR_CALL( CMRgraphFree(cmr, &dec->cograph) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographForest) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographCoforest) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographArcsReversed) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->reductions) );
  CMR_CALL( CMRchrmatFree(cmr, &dec->nestedMinorsMatrix) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsSequence) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsRowsOriginal) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsColumnsOriginal) );

  CMR_CALL( CMRfreeBlock(cmr, pdec) );

  return CMR_OKAY;
}

int CMRdecIsSum(CMR_DEC* dec, bool* plowerLeftNonzeros, bool* pupperRightNonzeros)
{
  assert(dec);

  if (dec->type == CMR_DEC_THREE_SUM)
  {
    if (plowerLeftNonzeros)
      assert("Not implemented for 3-sums." == 0);
    if (pupperRightNonzeros)
      assert("Not implemented for 3-sums." == 0);
    return 3;
  }
  else if (dec->type == CMR_DEC_TWO_SUM)
  {
    if (plowerLeftNonzeros)
      assert("Not implemented for 2-sums." == 0);
    if (pupperRightNonzeros)
      assert("Not implemented for 2-sums." == 0);
    return 2;
  }
  
  if (plowerLeftNonzeros)
    *plowerLeftNonzeros = 0;
  if (pupperRightNonzeros)
    *pupperRightNonzeros = 0;
  return dec->type == CMR_DEC_ONE_SUM ? 1 : 0;
}

bool CMRdecIsGraphicLeaf(CMR_DEC* dec)
{
  assert(dec);

  return dec->type == CMR_DEC_GRAPHIC;
}

bool CMRdecIsCographicLeaf(CMR_DEC* dec)
{
  assert(dec);

  return dec->type == CMR_DEC_COGRAPHIC;
}

CMR_DEC_TYPE CMRdecIsSpecialLeaf(CMR_DEC* dec, int* prepresentationMatrix)
{
  assert(dec);

  if (dec->type < CMR_DEC_SPECIAL_R10)
  {
    return 0;
  }
  else
  {
    if (prepresentationMatrix)
      *prepresentationMatrix = dec->flags & CMR_DEC_MASK_REPRESENTATION;
    return dec->type;
  }
}


bool CMRdecIsGraphic(CMR_DEC* dec)
{
  assert(dec);

  return dec->flags & CMR_DEC_IS_GRAPHIC;
}

bool CMRdecIsCographic(CMR_DEC* dec)
{
  assert(dec);
  
  return dec->flags & CMR_DEC_IS_COGRAPHIC;
}

bool CMRdecIsRegular(CMR_DEC* dec)
{
  assert(dec);

  return dec->flags & CMR_DEC_IS_REGULAR;
}

bool CMRdecNumRows(CMR_DEC* dec)
{
  assert(dec);

  return dec->numRows;
}

size_t* CMRdecRowsParent(CMR_DEC* dec)
{
  assert(dec);

  return dec->rowsParent;
}

bool CMRdecNumColumns(CMR_DEC* dec)
{
  assert(dec);

  return dec->numColumns;
}

size_t* CMRdecColumnsParent(CMR_DEC* dec)
{
  assert(dec);

  return dec->columnsParent;
}

CMR_ERROR CMRdecPrint(CMR* cmr, CMR_DEC* dec, FILE* stream, size_t indent, bool printMatrices, bool printGraphs,
  bool printReductions)
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
  
  fprintf(stream, "%ldx%ld ", dec->numRows, dec->numColumns);
  switch (dec->type)
  {
  case CMR_DEC_IRREGULAR:
    fprintf(stream, "irregular {");
  break;  
  case CMR_DEC_UNKNOWN:
    fprintf(stream, "unknown {");
  break;
  case CMR_DEC_ONE_SUM:
  case CMR_DEC_TWO_SUM:
  case CMR_DEC_THREE_SUM:
    fprintf(stream, "%d-sum with %ld children {", dec->type, dec->numChildren);
  break;
  case CMR_DEC_GRAPHIC:
    fprintf(stream, "graphic matrix with %d nodes and %d edges {", CMRgraphNumNodes(dec->graph), CMRgraphNumEdges(dec->graph));
  break;
  case CMR_DEC_COGRAPHIC:
    fprintf(stream, "cographic matrix with %d nodes and %d edges {", CMRgraphNumNodes(dec->cograph), CMRgraphNumEdges(dec->cograph));
  break;
  case CMR_DEC_PLANAR:
    assert(CMRgraphNumEdges(dec->graph) == CMRgraphNumEdges(dec->cograph));
    fprintf(stream, "planar matrix with %d nodes, %d faces and %d edges {", CMRgraphNumNodes(dec->graph),
      CMRgraphNumNodes(dec->cograph), CMRgraphNumEdges(dec->graph));
  break;
  case CMR_DEC_SERIES_PARALLEL:
    if (dec->numChildren)
      fprintf(stream, "matrix with %ld series-parallel reductions; 1 child {", dec->numReductions);
    else
      fprintf(stream, "series-parallel matrix {");
  break;
  case CMR_DEC_SPECIAL_R10:
    fprintf(stream, "matrix representing R10 {");
  break;
  case CMR_DEC_SPECIAL_FANO:
    fprintf(stream, "matrix representing F_7 {");
  break;
  case CMR_DEC_SPECIAL_FANO_DUAL:
    fprintf(stream, "matrix representing F_7^* {");
  break;
  case CMR_DEC_SPECIAL_K_5:
    fprintf(stream, "matrix representing K_5 {");
  break;
  case CMR_DEC_SPECIAL_K_5_DUAL:
    fprintf(stream, "matrix representing K_5^* {");
  break;
  case CMR_DEC_SPECIAL_K_3_3:
    fprintf(stream, "matrix representing K_{3,3} {");
  break;
  case CMR_DEC_SPECIAL_K_3_3_DUAL:
    fprintf(stream, "matrix representing K_{3,3}^* {");
  break;
  default:
    fprintf(stream, "invalid");
    fflush(stream);
    return CMR_ERROR_INVALID;
  break;
  }

  bool isFirst = true;
  if (dec->flags & CMR_DEC_IS_REGULAR)
  {
    fprintf(stream, "regular");
    isFirst = false;
  }
  if (dec->flags & CMR_DEC_IS_GRAPHIC)
  {
    fprintf(stream, "%sgraphic", isFirst ? "" : ",");
    isFirst = false;
  }
  if (dec->flags & CMR_DEC_IS_COGRAPHIC)
  {
    fprintf(stream, "%scographic", isFirst ? "" : ",");
    isFirst = false;
  }
  fprintf(stream, "}\n");

  if (printMatrices)
  {
    if (dec->matrix || dec->transpose)
    {
      for (int i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "Matrix:\n");
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
      for (int i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "No matrix.\n");
    }
  }
  if (printGraphs)
  {
    if (dec->graph)
    {
      for (int i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "Graph:\n");
      CMR_CALL( CMRgraphPrint(stream, dec->graph) );
    }
    if (dec->cograph)
    {
      for (int i = 0; i < indent; ++i)
        fputc(' ', stream);
      fprintf(stream, "Cograph:\n");
      CMR_CALL( CMRgraphPrint(stream, dec->cograph) );
    }
  }

  for (size_t c = 0; c < dec->numChildren; ++c)
    CMR_CALL( CMRdecPrint(cmr, dec->children[c], stream, indent + 2, printMatrices, printGraphs, printReductions) );

  return CMR_OKAY;
}


CMR_ERROR CMRdecCreate(CMR* cmr, CMR_DEC* parent, size_t numRows, size_t* rowsParent, size_t numColumns,
  size_t* columnsParent, CMR_DEC** pdec)
{
  assert(cmr);

  CMR_CALL( CMRallocBlock(cmr, pdec) );
  CMR_DEC* dec = *pdec;
  dec->type = CMR_DEC_UNKNOWN;
  dec->flags = 0;
  dec->matrix = NULL;
  dec->transpose = NULL;

  dec->parent = parent;
  dec->numChildren = 0;
  dec->children = NULL;

  dec->numRows = numRows;
  dec->rowsParent = NULL;
  if (rowsParent)
    CMR_CALL( CMRduplicateBlockArray(cmr, &dec->rowsParent, numRows, rowsParent) );
  dec->numColumns = numColumns;
  dec->columnsParent = NULL;
  if (columnsParent)
    CMR_CALL( CMRduplicateBlockArray(cmr, &dec->columnsParent, numColumns, columnsParent) );

  dec->graph = NULL;
  dec->graphForest = NULL;
  dec->graphCoforest = NULL;
  dec->graphArcsReversed = NULL;

  dec->cograph = NULL;
  dec->cographForest = NULL;
  dec->cographCoforest = NULL;
  dec->cographArcsReversed = NULL;

  dec->reductions = NULL;
  dec->numReductions = 0;

  dec->nestedMinorsMatrix = NULL;
  dec->nestedMinorsSequence = NULL;
  dec->nestedMinorsLength = 0;
  dec->nestedMinorsRowsOriginal = NULL;
  dec->nestedMinorsColumnsOriginal = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRdecInheritMatrices(CMR* cmr, CMR_DEC* node)
{
  assert(cmr);
  assert(node);

  /* First check if matrix or transpose is known. */
  if (node->matrix)
  {
    if (!node->transpose)
      CMR_CALL( CMRchrmatTranspose(cmr, node->matrix, &node->transpose) );
    return CMR_OKAY;
  }
  else if (node->transpose)
  {
    CMR_CALL( CMRchrmatTranspose(cmr, node->transpose, &node->matrix) );
    return CMR_OKAY;
  }

  assert(!node->matrix);
  assert(!node->transpose);
  assert(node->parent);

  /* We have to create the matrix. */
  CMR_CALL( CMRchrmatFilter(cmr, node->parent->matrix, node->numRows, node->rowsParent, node->numColumns,
    node->columnsParent, &node->matrix) );
  CMR_CALL( CMRchrmatTranspose(cmr, node->matrix, &node->transpose) );

  return CMR_OKAY;
}

CMR_ERROR CMRdecSetNumChildren(CMR* cmr, CMR_DEC* node, size_t numChildren)
{
  assert(cmr);
  assert(node);
  assert(numChildren >= 0);

  CMR_CALL( CMRreallocBlockArray(cmr, &node->children, numChildren) );
  for (size_t c = node->numChildren; c < numChildren; ++c)
    node->children[c] = NULL;
  node->numChildren = numChildren;

  return CMR_OKAY;
}

CMR_ERROR CMRdecComputeRegularity(CMR_DEC* node)
{
  assert(node);

  /* We mark it as graphic, cographic and regular and process the children first. */
  node->flags |= CMR_DEC_IS_GRAPHIC | CMR_DEC_IS_COGRAPHIC | CMR_DEC_IS_REGULAR;
  for (size_t c = 0; c < node->numChildren; ++c)
  {
    CMR_CALL( CMRdecComputeRegularity(node->children[c]) );

    CMR_DEC_FLAGS childFlags = node->children[c]->flags;
    if (!(childFlags & CMR_DEC_IS_GRAPHIC))
      node->flags &= ~CMR_DEC_IS_GRAPHIC;
    if (!(childFlags & CMR_DEC_IS_COGRAPHIC))
      node->flags &= ~CMR_DEC_IS_COGRAPHIC;
    if (!(childFlags & CMR_DEC_IS_REGULAR))
      node->flags &= ~CMR_DEC_IS_REGULAR;
  }

  /* We only need to set to false if appropriate. */
  switch (node->type)
  {
    case CMR_DEC_IRREGULAR:
    case CMR_DEC_UNKNOWN:
    case CMR_DEC_SPECIAL_FANO:
    case CMR_DEC_SPECIAL_FANO_DUAL:
      /* Set all to false. */
      node->flags &= ~(CMR_DEC_IS_GRAPHIC | CMR_DEC_IS_COGRAPHIC | CMR_DEC_IS_REGULAR);
    break;
    case CMR_DEC_ONE_SUM:
    case CMR_DEC_TWO_SUM:
    case CMR_DEC_THREE_SUM:
    case CMR_DEC_PLANAR:
      /* Do nothing: we simply inherit or have all flags true. */
    break;
    case CMR_DEC_GRAPHIC:
    case CMR_DEC_SPECIAL_K_5:
    case CMR_DEC_SPECIAL_K_3_3:
      node->flags &= ~(CMR_DEC_IS_COGRAPHIC);
    break;
    case CMR_DEC_COGRAPHIC:
    case CMR_DEC_SPECIAL_K_5_DUAL:
    case CMR_DEC_SPECIAL_K_3_3_DUAL:
      node->flags &= ~(CMR_DEC_IS_GRAPHIC);
    break;
    case CMR_DEC_SPECIAL_R10:
      node->flags &= ~(CMR_DEC_IS_GRAPHIC | CMR_DEC_IS_COGRAPHIC);
    break;
    case CMR_DEC_SERIES_PARALLEL:
    break;
    default:
      assert(false);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdecTranslateMinorToParent(CMR_DEC* node, CMR_MINOR* minor)
{
  assert(node);

  if (minor)
  {
    for (size_t pivot = 0; pivot < minor->numPivots; ++pivot)
    {
      minor->pivotRows[pivot] = node->rowsParent[minor->pivotRows[pivot]];
      minor->pivotColumns[pivot] = node->columnsParent[minor->pivotColumns[pivot]];
    }

    assert(minor->remainingSubmatrix);
    for (size_t row = 0; row < minor->remainingSubmatrix->numRows; ++row)
      minor->remainingSubmatrix->rows[row] = node->rowsParent[minor->remainingSubmatrix->rows[row]];
    for (size_t column = 0; column < minor->remainingSubmatrix->numColumns; ++column)
      minor->remainingSubmatrix->columns[column] = node->columnsParent[minor->remainingSubmatrix->columns[column]];
  }

  return CMR_OKAY;
}

CMR_ERROR CMRdecApplySeparation(CMR* cmr, CMR_DEC* dec, CMR_SEPA* sepa)
{
  assert(cmr);
  assert(dec);
  assert(sepa);

  if (CMRsepaRank(sepa) == 1)
  {
    dec->type = CMR_DEC_TWO_SUM;
    CMR_CALL( CMRdecSetNumChildren(cmr, dec, 2) );
    unsigned char rankBottomLeft = CMRsepaRankBottomLeft(sepa);
    unsigned char rankTopRight = CMRsepaRankTopRight(sepa);

    CMRdbgMsg(4, "Ranks are %d and %d.\n", rankBottomLeft, rankTopRight);

    for (size_t child = 0; child < 2; ++child)
    {
      unsigned char numExtraRows = child == 0 ? rankBottomLeft : rankTopRight;
      unsigned char numExtraColumns = child == 0 ? rankTopRight : rankBottomLeft;

      CMRdbgMsg(4, "Child %d has %d+%d rows and %d+%d columns.\n", child, sepa->numRows[child], numExtraRows,
        sepa->numColumns[child], numExtraColumns);

      size_t* extraRows = child == 0 ? sepa->extraRows0 : sepa->extraRows1;
      size_t* extraColumns = child == 0 ? sepa->extraColumns0 : sepa->extraColumns1;

      CMR_CALL( CMRdecCreate(cmr, dec, sepa->numRows[child] + numExtraRows, NULL,
        sepa->numColumns[child] + numExtraColumns, NULL, &dec->children[child]) );

      CMR_CALL( CMRallocBlockArray(cmr, &dec->children[child]->rowsParent, dec->children[child]->numRows) );
      for (size_t row = 0; row < sepa->numRows[child]; ++row)
        dec->children[child]->rowsParent[row] = sepa->rows[child][row];
      for (unsigned char extra = 0; extra < numExtraRows; ++extra)
        dec->children[child]->rowsParent[sepa->numRows[child] + extra] = extraRows[extra];

      CMR_CALL( CMRallocBlockArray(cmr, &dec->children[child]->columnsParent, dec->children[child]->numColumns) );
      for (size_t column = 0; column < sepa->numColumns[child]; ++column)
        dec->children[child]->columnsParent[column] = sepa->columns[child][column];
      for (unsigned char extra = 0; extra < numExtraColumns; ++extra)
        dec->children[child]->columnsParent[sepa->numColumns[child] + extra] = extraColumns[extra];
      
      CMR_CALL( CMRdecInheritMatrices(cmr, dec->children[child]) );
    }
  }
  else
  {
    assert("CMRdecApplySeparation only implemented for 2-sums." == 0);
  }

  return CMR_OKAY;
}

#include "decomposition_internal.h"

#include "matrix_internal.h"

CMR_ERROR CMRtudecFree(CMR* cmr, CMR_TU_DEC** pdec)
{
  assert(cmr);
  assert(pdec);

  CMR_TU_DEC* dec = *pdec;
  if (!dec)
    return CMR_OKAY;

  for (size_t c = 0; c < dec->numChildren; ++c)
    CMR_CALL( CMRtudecFree(cmr, &dec->children[c]) );

  CMR_CALL( CMRchrmatFree(cmr, &dec->matrix) );
  CMR_CALL( CMRchrmatFree(cmr, &dec->transpose) );
  CMR_CALL( CMRfreeBlockArray(cmr, dec->rowsParent) );
  CMR_CALL( CMRfreeBlockArray(cmr, dec->columnsParent) );
  CMR_CALL( CMRfreeBlockArray(cmr, dec->rowElements) );
  CMR_CALL( CMRfreeBlockArray(cmr, dec->columnElements) );
  CMR_CALL( CMRgraphFree(cmr, &dec->graph) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->edgeElements) );
  CMR_CALL( CMRgraphFree(cmr, &dec->cograph) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->coedgeElements) );

  CMR_CALL( CMRfreeBlock(cmr, pdec) );

  return CMR_OKAY;
}

int CMRtudecIsSum(CMR_TU_DEC* dec, bool* plowerLeftNonzeros, bool* pupperRightNonzeros)
{
  assert(dec);

  if (dec->type == CMR_TU_DEC_THREE_SUM)
  {
    if (plowerLeftNonzeros)
      assert(false);
    if (pupperRightNonzeros)
      assert(false);
    return 3;
  }
  else if (dec->type == CMR_TU_DEC_TWO_SUM)
  {
    if (plowerLeftNonzeros)
      assert(false);
    if (pupperRightNonzeros)
      assert(false);
    return 2;
  }
  
  if (plowerLeftNonzeros)
    *plowerLeftNonzeros = 0;
  if (pupperRightNonzeros)
    *pupperRightNonzeros = 0;
  return dec->type == CMR_TU_DEC_ONE_SUM ? 1 : 0;
}

bool CMRtudecIsPivotSequence(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->type == CMR_TU_DEC_PIVOTS;
}

bool CMRtudecIsSubmatrix(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->type == CMR_TU_DEC_SUBMATRIX;
}

bool CMRtudecIsGraphicLeaf(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->type == CMR_TU_DEC_GRAPHIC;
}

bool CMRtudecIsCographicLeaf(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->type == CMR_TU_DEC_COGRAPHIC;
}

CMR_TU_DEC_TYPE CMRtudecIsSpecialLeaf(CMR_TU_DEC* dec, int* prepresentationMatrix)
{
  assert(dec);

  if (dec->type < CMR_TU_DEC_SPECIAL_R10)
  {
    return 0;
  }
  else
  {
    if (prepresentationMatrix)
      *prepresentationMatrix = dec->flags & CMR_TU_DEC_MASK_REPRESENTATION;
    return dec->type;
  }
}


bool CMRtudecIsGraphic(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->flags & CMR_TU_DEC_IS_GRAPHIC;
}

bool CMRtudecIsCographic(CMR_TU_DEC* dec)
{
  assert(dec);
  
  return dec->flags & CMR_TU_DEC_IS_COGRAPHIC;
}

bool CMRtudecIsRegular(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->flags & CMR_TU_DEC_IS_REGULAR;
}

bool CMRtudecNumRows(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->numRows;
}

CMR_ELEMENT* CMRtudecRowElements(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->rowElements;
}

size_t* CMRtudecRowsParent(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->rowsParent;
}

bool CMRtudecNumColumns(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->numColumns;
}

CMR_ELEMENT* CMRtudecColumnElements(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->columnElements;
}

size_t* CMRtudecColumnsParent(CMR_TU_DEC* dec)
{
  assert(dec);

  return dec->columnsParent;
}

CMR_ERROR CMRtudecPrint(FILE* stream, CMR_TU_DEC* dec, size_t indent)
{
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
  case CMR_TU_DEC_IRREGULAR:
    fprintf(stream, "irregular {");
  break;  
  case CMR_TU_DEC_UNKNOWN:
    fprintf(stream, "unknown {");
  break;
  case CMR_TU_DEC_ONE_SUM:
  case CMR_TU_DEC_TWO_SUM:
  case CMR_TU_DEC_THREE_SUM:
    fprintf(stream, "%d-sum with %ld children {", dec->type, dec->numChildren);
  break;
  case CMR_TU_DEC_SUBMATRIX:
    fprintf(stream, "submatrix {");
  break;
  case CMR_TU_DEC_PIVOTS:
    fprintf(stream, "%ld pivots {", dec->numPivots);
  break;
  case CMR_TU_DEC_GRAPHIC:
    fprintf(stream, "graphic with %d nodes and %d edges {", CMRgraphNumNodes(dec->graph), CMRgraphNumEdges(dec->graph));
  break;
  case CMR_TU_DEC_COGRAPHIC:
    fprintf(stream, "cographic with %d nodes and %d edges {", CMRgraphNumNodes(dec->graph), CMRgraphNumEdges(dec->graph));
  break;
  case CMR_TU_DEC_PLANAR:
    assert(CMRgraphNumEdges(dec->graph) == CMRgraphNumEdges(dec->cograph));
    fprintf(stream, "planar with %d nodes, %d faces and %d edges {", CMRgraphNumNodes(dec->graph),
      CMRgraphNumNodes(dec->cograph), CMRgraphNumEdges(dec->graph));
  break;
  case CMR_TU_DEC_SPECIAL_R10:
    fprintf(stream, "R10 {");
  break;
  case CMR_TU_DEC_SPECIAL_FANO:
    fprintf(stream, "F_7 {");
  break;
  case CMR_TU_DEC_SPECIAL_FANO_DUAL:
    fprintf(stream, "F_7^* {");
  break;
  case CMR_TU_DEC_SPECIAL_K_5:
    fprintf(stream, "K_5 {");
  break;
  case CMR_TU_DEC_SPECIAL_K_5_DUAL:
    fprintf(stream, "K_5^* {");
  break;
  case CMR_TU_DEC_SPECIAL_K_3_3:
    fprintf(stream, "K_{3,3} {");
  break;
  case CMR_TU_DEC_SPECIAL_K_3_3_DUAL:
    fprintf(stream, "K_{3,3}^* {");
  break;
  default:
    fprintf(stream, "invalid");
    fflush(stream);
    return CMR_ERROR_INVALID;
  break;
  }

  bool isFirst = true;
  if (dec->flags & CMR_TU_DEC_IS_REGULAR)
  {
    fprintf(stream, "regular");
    isFirst = false;
  }
  if (dec->flags & CMR_TU_DEC_IS_GRAPHIC)
  {
    fprintf(stream, "%sgraphic", isFirst ? "" : ",");
    isFirst = false;
  }
  if (dec->flags & CMR_TU_DEC_IS_COGRAPHIC)
  {
    fprintf(stream, "%scographic", isFirst ? "" : ",");
    isFirst = false;
  }
  fprintf(stream, "}\n");

  for (size_t c = 0; c < dec->numChildren; ++c)
    CMR_CALL( CMRtudecPrint(stream, dec->children[c], indent + 2) );

  return CMR_OKAY;
}


CMR_ERROR CMRtudecCreate(CMR* cmr, CMR_TU_DEC* parent, size_t numRows, size_t* rowsParent, size_t numColumns,
  size_t* columnsParent, CMR_TU_DEC** pdec)
{
  assert(cmr);

  CMR_CALL( CMRallocBlock(cmr, pdec) );
  CMR_TU_DEC* dec = *pdec;
  dec->type = CMR_TU_DEC_UNKNOWN;
  dec->flags = 0;
  dec->matrix = NULL;
  dec->transpose = NULL;

  dec->numRows = numRows;
  dec->rowElements = NULL;
  dec->rowsParent = NULL;
  if (rowsParent)
    CMR_CALL( CMRduplicateBlockArray(cmr, &dec->rowsParent, numRows, rowsParent) );
  dec->numColumns = numColumns;
  dec->columnElements = NULL;
  dec->columnsParent = NULL;
  if (columnsParent)
    CMR_CALL( CMRduplicateBlockArray(cmr, &dec->columnsParent, numColumns, columnsParent) );

  dec->numPivots = 0;
  dec->graph = NULL;
  dec->edgeElements = NULL;
  dec->cograph = NULL;
  dec->coedgeElements = NULL;
  dec->parent = parent;
  dec->numChildren = 0;
  dec->children = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRtudecInheritElements(CMR* cmr, CMR_TU_DEC* node)
{
  assert(cmr);
  assert(node);
  assert(node->parent);

  /* We assume that there is no row/column element information, yet. */
  assert(!node->rowElements);
  assert(!node->columnElements);

  CMR_CALL( CMRallocBlockArray(cmr, &node->rowElements, node->numRows) );
  for (size_t row = 0; row < node->numRows; ++row)
    node->rowElements[row] = node->parent->rowElements[node->rowsParent[row]];

  CMR_CALL( CMRallocBlockArray(cmr, &node->columnElements, node->numColumns) );
  for (size_t column = 0; column < node->numColumns; ++column)
    node->columnElements[column] = node->parent->columnElements[node->columnsParent[column]];

  return CMR_OKAY;
}

CMR_ERROR CMRtudecInheritMatrices(CMR* cmr, CMR_TU_DEC* node)
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

CMR_ERROR CMRtudecSetNumChildren(CMR* cmr, CMR_TU_DEC* node, size_t numChildren)
{
  assert(cmr);
  assert(node);
  assert(numChildren >= 0);

  node->numChildren = numChildren;
  CMR_CALL( CMRreallocBlockArray(cmr, &node->children, numChildren) );
  for (size_t c = 0; c < numChildren; ++c)
    node->children[c] = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRtudecComputeRegularity(CMR_TU_DEC* node)
{
  assert(node);

  /* We mark it as graphic, cographic and regular and process the children first. */
  node->flags |= CMR_TU_DEC_IS_GRAPHIC | CMR_TU_DEC_IS_COGRAPHIC | CMR_TU_DEC_IS_REGULAR;
  for (size_t c = 0; c < node->numChildren; ++c)
  {
    CMR_CALL( CMRtudecComputeRegularity(node->children[c]) );

    CMR_TU_DEC_FLAGS childFlags = node->children[c]->flags;
    if (!(childFlags & CMR_TU_DEC_IS_GRAPHIC))
      node->flags &= ~CMR_TU_DEC_IS_GRAPHIC;
    if (!(childFlags & CMR_TU_DEC_IS_COGRAPHIC))
      node->flags &= ~CMR_TU_DEC_IS_COGRAPHIC;
    if (!(childFlags & CMR_TU_DEC_IS_REGULAR))
      node->flags &= ~CMR_TU_DEC_IS_REGULAR;
  }

  /* We only need to set to false if appropriate. */
  switch (node->type)
  {
    case CMR_TU_DEC_IRREGULAR:
    case CMR_TU_DEC_UNKNOWN:
    case CMR_TU_DEC_SUBMATRIX:
    case CMR_TU_DEC_SPECIAL_FANO:
    case CMR_TU_DEC_SPECIAL_FANO_DUAL:
      /* Set all to false. */
      node->flags &= ~(CMR_TU_DEC_IS_GRAPHIC | CMR_TU_DEC_IS_COGRAPHIC | CMR_TU_DEC_IS_REGULAR);
    break;
    case CMR_TU_DEC_ONE_SUM:
    case CMR_TU_DEC_TWO_SUM:
    case CMR_TU_DEC_THREE_SUM:
    case CMR_TU_DEC_PIVOTS:
    case CMR_TU_DEC_PLANAR:
      /* Do nothing: we simply inherit or have all flags true. */
    break;
    case CMR_TU_DEC_GRAPHIC:
    case CMR_TU_DEC_SPECIAL_K_5:
    case CMR_TU_DEC_SPECIAL_K_3_3:
      node->flags &= ~(CMR_TU_DEC_IS_COGRAPHIC);
    break;
    case CMR_TU_DEC_COGRAPHIC:
    case CMR_TU_DEC_SPECIAL_K_5_DUAL:
    case CMR_TU_DEC_SPECIAL_K_3_3_DUAL:
      node->flags &= ~(CMR_TU_DEC_IS_GRAPHIC);
    break;
    case CMR_TU_DEC_SPECIAL_R10:
      node->flags &= ~(CMR_TU_DEC_IS_GRAPHIC | CMR_TU_DEC_IS_COGRAPHIC);
    break;
    default:
      assert(false);
  }

  return CMR_OKAY;
}


#include "decomposition_internal.h"

#include "matrix_internal.h"

TU_ERROR TUdecFree(TU* tu, TU_DEC** pdec)
{
  assert(tu);
  assert(pdec);

  TU_DEC* dec = *pdec;
  if (!dec)
    return TU_OKAY;

  for (size_t c = 0; c < dec->numChildren; ++c)
    TU_CALL( TUdecFree(tu, &dec->children[c]) );

  TU_CALL( TUchrmatFree(tu, &dec->matrix) );
  TU_CALL( TUchrmatFree(tu, &dec->transpose) );
  TU_CALL( TUfreeBlockArray(tu, dec->rowsParent) );
  TU_CALL( TUfreeBlockArray(tu, dec->columnsParent) );
  TU_CALL( TUfreeBlockArray(tu, dec->rowElements) );
  TU_CALL( TUfreeBlockArray(tu, dec->columnElements) );
  TU_CALL( TUgraphFree(tu, &dec->graph) );
  TU_CALL( TUfreeBlockArray(tu, &dec->edgeElements) );
  TU_CALL( TUgraphFree(tu, &dec->cograph) );
  TU_CALL( TUfreeBlockArray(tu, &dec->coedgeElements) );

  TU_CALL( TUfreeBlock(tu, pdec) );

  return TU_OKAY;
}

int TUdecIsSum(TU_DEC* dec, bool* plowerLeftNonzeros, bool* pupperRightNonzeros)
{
  assert(dec);

  if (dec->type == TU_DEC_THREE_SUM)
  {
    if (plowerLeftNonzeros)
      assert(false);
    if (pupperRightNonzeros)
      assert(false);
    return 3;
  }
  else if (dec->type == TU_DEC_TWO_SUM)
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
  return dec->type == TU_DEC_ONE_SUM ? 1 : 0;
}

bool TUdecIsPivotSequence(TU_DEC* dec)
{
  assert(dec);

  return dec->type == TU_DEC_PIVOTS;
}

bool TUdecIsSubmatrix(TU_DEC* dec)
{
  assert(dec);

  return dec->type == TU_DEC_SUBMATRIX;
}

bool TUdecIsGraphicLeaf(TU_DEC* dec)
{
  assert(dec);

  return dec->type == TU_DEC_GRAPHIC;
}

bool TUdecIsCographicLeaf(TU_DEC* dec)
{
  assert(dec);

  return dec->type == TU_DEC_COGRAPHIC;
}

TU_DEC_TYPE TUdecIsSpecialLeaf(TU_DEC* dec, int* prepresentationMatrix)
{
  assert(dec);

  if (dec->type < TU_DEC_SPECIAL_R10)
  {
    return 0;
  }
  else
  {
    if (prepresentationMatrix)
      *prepresentationMatrix = dec->flags & TU_DEC_MASK_REPRESENTATION;
    return dec->type;
  }
}


bool TUdecIsGraphic(TU_DEC* dec)
{
  assert(dec);

  return dec->flags & TU_DEC_IS_GRAPHIC;
}

bool TUdecIsCographic(TU_DEC* dec)
{
  assert(dec);
  
  return dec->flags & TU_DEC_IS_COGRAPHIC;
}

bool TUdecIsRegular(TU_DEC* dec)
{
  assert(dec);

  return dec->flags & TU_DEC_IS_REGULAR;
}

bool TUdecNumRows(TU_DEC* dec)
{
  assert(dec);

  return dec->numRows;
}

TU_ELEMENT* TUdecRowElements(TU_DEC* dec)
{
  assert(dec);

  return dec->rowElements;
}

size_t* TUdecRowsParent(TU_DEC* dec)
{
  assert(dec);

  return dec->rowsParent;
}

bool TUdecNumColumns(TU_DEC* dec)
{
  assert(dec);

  return dec->numColumns;
}

TU_ELEMENT* TUdecColumnElements(TU_DEC* dec)
{
  assert(dec);

  return dec->columnElements;
}

size_t* TUdecColumnsParent(TU_DEC* dec)
{
  assert(dec);

  return dec->columnsParent;
}

TU_ERROR TUdecPrint(FILE* stream, TU_DEC* dec, size_t indent)
{
  assert(stream);

  /* Indent. */
  for (size_t i = 0; i < indent; ++i)
    fputc(' ', stream);

  if (!dec)
  {
    fprintf(stream, "<NULL>\n");
    return TU_OKAY;
  }
  
  fprintf(stream, "%ldx%ld ", dec->numRows, dec->numColumns);
  switch (dec->type)
  {
  case TU_DEC_IRREGULAR:
    fprintf(stream, "irregular {");
  break;  
  case TU_DEC_UNKNOWN:
    fprintf(stream, "unknown {");
  break;
  case TU_DEC_ONE_SUM:
  case TU_DEC_TWO_SUM:
  case TU_DEC_THREE_SUM:
    fprintf(stream, "%d-sum with %ld children {", dec->type, dec->numChildren);
  break;
  case TU_DEC_SUBMATRIX:
    fprintf(stream, "submatrix {");
  break;
  case TU_DEC_PIVOTS:
    fprintf(stream, "%ld pivots {", dec->numPivots);
  break;
  case TU_DEC_GRAPHIC:
    fprintf(stream, "graphic with %d nodes and %d edges {", TUgraphNumNodes(dec->graph), TUgraphNumEdges(dec->graph));
  break;
  case TU_DEC_COGRAPHIC:
    fprintf(stream, "cographic with %d nodes and %d edges {", TUgraphNumNodes(dec->graph), TUgraphNumEdges(dec->graph));
  break;
  case TU_DEC_PLANAR:
    assert(TUgraphNumEdges(dec->graph) == TUgraphNumEdges(dec->cograph));
    fprintf(stream, "planar with %d nodes, %d faces and %d edges {", TUgraphNumNodes(dec->graph),
      TUgraphNumNodes(dec->cograph), TUgraphNumEdges(dec->graph));
  break;
  case TU_DEC_SPECIAL_R10:
    fprintf(stream, "R10 {");
  break;
  case TU_DEC_SPECIAL_FANO:
    fprintf(stream, "F_7 {");
  break;
  case TU_DEC_SPECIAL_FANO_DUAL:
    fprintf(stream, "F_7^* {");
  break;
  case TU_DEC_SPECIAL_K_5:
    fprintf(stream, "K_5 {");
  break;
  case TU_DEC_SPECIAL_K_5_DUAL:
    fprintf(stream, "K_5^* {");
  break;
  case TU_DEC_SPECIAL_K_3_3:
    fprintf(stream, "K_{3,3} {");
  break;
  case TU_DEC_SPECIAL_K_3_3_DUAL:
    fprintf(stream, "K_{3,3}^* {");
  break;
  default:
    fprintf(stream, "invalid");
    fflush(stream);
    return TU_ERROR_INVALID;
  break;
  }

  bool isFirst = true;
  if (dec->flags & TU_DEC_IS_REGULAR)
  {
    fprintf(stream, "regular");
    isFirst = false;
  }
  if (dec->flags & TU_DEC_IS_GRAPHIC)
  {
    fprintf(stream, "%sgraphic", isFirst ? "" : ",");
    isFirst = false;
  }
  if (dec->flags & TU_DEC_IS_COGRAPHIC)
  {
    fprintf(stream, "%scographic", isFirst ? "" : ",");
    isFirst = false;
  }
  fprintf(stream, "}\n");

  for (size_t c = 0; c < dec->numChildren; ++c)
    TU_CALL( TUdecPrint(stream, dec->children[c], indent + 2) );

  return TU_OKAY;
}


TU_ERROR TUdecCreate(TU* tu, TU_DEC* parent, size_t numRows, size_t* rowsParent, size_t numColumns,
  size_t* columnsParent, TU_DEC** pdec)
{
  assert(tu);

  TU_CALL( TUallocBlock(tu, pdec) );
  TU_DEC* dec = *pdec;
  dec->type = TU_DEC_UNKNOWN;
  dec->flags = 0;
  dec->matrix = NULL;
  dec->transpose = NULL;

  dec->numRows = numRows;
  dec->rowElements = NULL;
  dec->rowsParent = NULL;
  if (rowsParent)
    TU_CALL( TUduplicateBlockArray(tu, &dec->rowsParent, numRows, rowsParent) );
  dec->numColumns = numColumns;
  dec->columnElements = NULL;
  dec->columnsParent = NULL;
  if (columnsParent)
    TU_CALL( TUduplicateBlockArray(tu, &dec->columnsParent, numColumns, columnsParent) );

  dec->numPivots = 0;
  dec->graph = NULL;
  dec->edgeElements = NULL;
  dec->cograph = NULL;
  dec->coedgeElements = NULL;
  dec->parent = parent;
  dec->numChildren = 0;
  dec->children = NULL;

  return TU_OKAY;
}

TU_ERROR TUdecInheritElements(TU* tu, TU_DEC* node)
{
  assert(tu);
  assert(node);
  assert(node->parent);

  /* We assume that there is no row/column element information, yet. */
  assert(!node->rowElements);
  assert(!node->columnElements);

  TU_CALL( TUallocBlockArray(tu, &node->rowElements, node->numRows) );
  for (size_t row = 0; row < node->numRows; ++row)
    node->rowElements[row] = node->parent->rowElements[node->rowsParent[row]];

  TU_CALL( TUallocBlockArray(tu, &node->columnElements, node->numColumns) );
  for (size_t column = 0; column < node->numColumns; ++column)
    node->columnElements[column] = node->parent->columnElements[node->columnsParent[column]];

  return TU_OKAY;
}

TU_ERROR TUdecInheritMatrices(TU* tu, TU_DEC* node)
{
  assert(tu);
  assert(node);

  /* First check if matrix or transpose is known. */
  if (node->matrix)
  {
    if (!node->transpose)
      TU_CALL( TUchrmatTranspose(tu, node->matrix, &node->transpose) );
    return TU_OKAY;
  }
  else if (node->transpose)
  {
    TU_CALL( TUchrmatTranspose(tu, node->transpose, &node->matrix) );
    return TU_OKAY;
  }

  assert(!node->matrix);
  assert(!node->transpose);
  assert(node->parent);

  /* We have to create the matrix. */
  TU_CALL( TUchrmatFilter(tu, node->parent->matrix, node->numRows, node->rowsParent, node->numColumns,
    node->columnsParent, &node->matrix) );
  TU_CALL( TUchrmatTranspose(tu, node->matrix, &node->transpose) );

  return TU_OKAY;
}

TU_ERROR TUdecSetNumChildren(TU* tu, TU_DEC* node, size_t numChildren)
{
  assert(tu);
  assert(node);
  assert(numChildren >= 0);

  node->numChildren = numChildren;
  TU_CALL( TUreallocBlockArray(tu, &node->children, numChildren) );
  for (size_t c = 0; c < numChildren; ++c)
    node->children[c] = NULL;

  return TU_OKAY;
}

TU_ERROR TUdecComputeRegularity(TU_DEC* node)
{
  assert(node);

  /* We mark it as graphic, cographic and regular and process the children first. */
  node->flags |= TU_DEC_IS_GRAPHIC | TU_DEC_IS_COGRAPHIC | TU_DEC_IS_REGULAR;
  for (size_t c = 0; c < node->numChildren; ++c)
  {
    TU_CALL( TUdecComputeRegularity(node->children[c]) );

    TU_DEC_FLAGS childFlags = node->children[c]->flags;
    if (!(childFlags & TU_DEC_IS_GRAPHIC))
      node->flags &= ~TU_DEC_IS_GRAPHIC;
    if (!(childFlags & TU_DEC_IS_COGRAPHIC))
      node->flags &= ~TU_DEC_IS_COGRAPHIC;
    if (!(childFlags & TU_DEC_IS_REGULAR))
      node->flags &= ~TU_DEC_IS_REGULAR;
  }

  /* We only need to set to false if appropriate. */
  switch (node->type)
  {
    case TU_DEC_IRREGULAR:
    case TU_DEC_UNKNOWN:
    case TU_DEC_SUBMATRIX:
    case TU_DEC_SPECIAL_FANO:
    case TU_DEC_SPECIAL_FANO_DUAL:
      /* Set all to false. */
      node->flags &= ~(TU_DEC_IS_GRAPHIC | TU_DEC_IS_COGRAPHIC | TU_DEC_IS_REGULAR);
    break;
    case TU_DEC_ONE_SUM:
    case TU_DEC_TWO_SUM:
    case TU_DEC_THREE_SUM:
    case TU_DEC_PIVOTS:
    case TU_DEC_PLANAR:
      /* Do nothing: we simply inherit or have all flags true. */
    break;
    case TU_DEC_GRAPHIC:
    case TU_DEC_SPECIAL_K_5:
    case TU_DEC_SPECIAL_K_3_3:
      node->flags &= ~(TU_DEC_IS_COGRAPHIC);
    break;
    case TU_DEC_COGRAPHIC:
    case TU_DEC_SPECIAL_K_5_DUAL:
    case TU_DEC_SPECIAL_K_3_3_DUAL:
      node->flags &= ~(TU_DEC_IS_GRAPHIC);
    break;
    case TU_DEC_SPECIAL_R10:
      node->flags &= ~(TU_DEC_IS_GRAPHIC | TU_DEC_IS_COGRAPHIC);
    break;
    default:
      assert(false);
  }

  return TU_OKAY;
}


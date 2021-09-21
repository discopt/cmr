// #define CMR_DEBUG /* Uncomment to debug this file. */

#include "env_internal.h"
#include "dec_internal.h"
#include "matrix_internal.h"

#include <string.h>

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
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsSequenceNumColumns) );
  CMR_CALL( CMRfreeBlockArray(cmr, &dec->nestedMinorsSequenceNumRows) );
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
  dec->nestedMinorsSequenceNumRows= NULL;
  dec->nestedMinorsSequenceNumColumns = NULL;
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

  unsigned char rankBottomLeft = CMRsepaRankBottomLeft(sepa);
  unsigned char rankTopRight = CMRsepaRankTopRight(sepa);
  CMRdbgMsg(4, "Ranks are %d and %d.\n", rankBottomLeft, rankTopRight);

  if (CMRsepaRank(sepa) == 1)
  {
    dec->type = CMR_DEC_TWO_SUM;
    CMR_CALL( CMRdecSetNumChildren(cmr, dec, 2) );
    

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

CMR_ERROR CMRdecPrintSequenceNested3ConnectedMinors(CMR* cmr, CMR_DEC* dec, FILE* stream)
{
  assert(cmr);
  assert(dec);
  assert(stream);

  if (!dec->nestedMinorsMatrix)
    return CMR_OKAY;

  bool isComplete = dec->nestedMinorsSequenceNumColumns[dec->nestedMinorsLength-1] == dec->matrix->numColumns
    && dec->nestedMinorsSequenceNumRows[dec->nestedMinorsLength-1] == dec->matrix->numRows;

  for (int i = 0; i < 3; ++i)
  {
    fprintf(stream, "%*s", 3+ 2*i, "");
    for (size_t column = i; column < dec->matrix->numColumns; column += 3)
    {
      char buffer[16];
      strcpy(buffer, CMRelementString(dec->nestedMinorsColumnsOriginal[column], 0));
      size_t length = strlen(buffer);
      int leftPadding = (6-length) / 2;
      int rightPadding = 6-length - leftPadding;
      fprintf(stream, "%*s%s%*s", leftPadding, "", buffer, rightPadding, "");
    }
    fputs("\n", stream);
  }

  fputs("    ", stream);
  for (size_t column = 0; column < dec->matrix->numColumns; ++column)
    fputs(column == 0 ? "+-" : "--", stream);
  fprintf(stream, "%s\n", isComplete ? "+" : "");

  size_t maxNested = 0;
  for (size_t row = 0; row < dec->matrix->numRows; ++row)
  {
    /* maxNested is the first minor with #rows > current row. */
    bool increased = false;
    while (maxNested < dec->nestedMinorsLength && dec->nestedMinorsSequenceNumRows[maxNested] <= row)
    {
      increased = true;
      ++maxNested;
    }

    if (increased)
    {
      fprintf(stream, "    ");
      for (size_t column = 0; column < dec->matrix->numColumns; ++column)
      {
        char separator = ' ';
        for (size_t j = maxNested; j < dec->nestedMinorsLength; ++j)
        {
          if (dec->nestedMinorsSequenceNumColumns[j] == column)
            separator = '|';
        }

        if (column < dec->nestedMinorsSequenceNumColumns[maxNested-1])
          fprintf(stream, column == 0 ? "+-" : "--");
        else if (column == dec->nestedMinorsSequenceNumColumns[maxNested-1])
          fprintf(stream, "+ ");
        else
          fprintf(stream, "%c ", separator);
      }
      fprintf(stream, "%s\n", isComplete ? "|" : "");
    }
    
    /* Row label */
    char buffer[16];
    strcpy(buffer, CMRelementString(dec->nestedMinorsRowsOriginal[row], 0));
    size_t length = strlen(buffer);
    int rightPadding = 4-length;
    fprintf(stream, "%s%*s", buffer, rightPadding, "");

    size_t first = dec->nestedMinorsMatrix->rowSlice[row];
    size_t beyond = dec->nestedMinorsMatrix->rowSlice[row + 1];
    size_t e = first;
    for (size_t column = 0; column < dec->matrix->numColumns; ++column)
    {
      char value;
      if (e < beyond && dec->nestedMinorsMatrix->entryColumns[e] == column)
      {
        value = '1';
        ++e;
      }
      else
        value = '0';
      char separator = column > 0 ? ' ' : '|';
      for (size_t j = maxNested; j < dec->nestedMinorsLength; ++j)
        if (dec->nestedMinorsSequenceNumColumns[j] == column)
          separator = '|';
      fprintf(stream, "%c%c", separator, value);
    }
    fprintf(stream, "%s\n", isComplete ? "|" : "");
  }

  if (isComplete)
  {
    fputs("    ", stream);
    for (size_t column = 0; column < dec->matrix->numColumns; ++column)
      fputs(column == 0 ? "+-" : "--", stream);
    fputs("+\n", stream);
  }

  return CMR_OKAY;
}


static
char* consistencyNested(
  CMR_DEC* dec  /**< Decomposition. */
)
{
  assert(dec);

  if (!dec->nestedMinorsMatrix)
    return NULL;

  if (!dec->matrix)
    return CMRconsistencyMessage("nested minor matrix exists, but matrix is missing.");
  if (!dec->nestedMinorsRowsOriginal)
    return CMRconsistencyMessage("nested minor matrix exists, row mapping is missing.");
  if (!dec->nestedMinorsColumnsOriginal)
    return CMRconsistencyMessage("nested minor matrix exists, column mapping is missing.");
  if (dec->nestedMinorsMatrix->numRows != dec->matrix->numRows)
  {
    return CMRconsistencyMessage("nested minor matrix has %ld rows, but matrix has %ld rows.",
      dec->nestedMinorsMatrix->numRows, dec->matrix->numRows);
  }
  if (dec->nestedMinorsMatrix->numColumns != dec->matrix->numColumns)
  {
    return CMRconsistencyMessage("nested minor matrix has %ld columns, but matrix has %ld columns.",
      dec->nestedMinorsMatrix->numColumns, dec->matrix->numColumns);
  }

  size_t numRows = dec->matrix->numRows;
  size_t numColumns = dec->matrix->numColumns;
  
  char** dense = (char**) malloc( numRows * sizeof(char*) );
  for (size_t row = 0; row < numRows; ++row)
  {
    dense[row] = (char*) malloc( numColumns * sizeof(char) );
    for (size_t column = 0; column < numColumns; ++column)
      dense[row][column] = 0;
    size_t first = dec->nestedMinorsMatrix->rowSlice[row];
    size_t beyond = dec->nestedMinorsMatrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
      dense[row][dec->nestedMinorsMatrix->entryColumns[e]] = dec->nestedMinorsMatrix->entryValues[e];
  }

  CMR_ELEMENT* rowElements = (CMR_ELEMENT*) malloc( numRows * sizeof(CMR_ELEMENT) );
  for (size_t row = 0; row < numRows; ++row)
    rowElements[row] = dec->nestedMinorsRowsOriginal[row];
  CMR_ELEMENT* columnElements = (CMR_ELEMENT*) malloc( numColumns * sizeof(CMR_ELEMENT) );
  for (size_t column = 0; column < numColumns; ++column)
    columnElements[column] = dec->nestedMinorsColumnsOriginal[column];

  while (true)
  {
    /* Find row that is labeled as a column. */
    size_t pivotRow = SIZE_MAX;
    for (size_t row = 0; row < numRows; ++row)
    {
      if (CMRelementIsColumn(rowElements[row]))
      {
        pivotRow = row;
        break;
      }
    }

    if (pivotRow == SIZE_MAX)
      break;

    /* In that row, find a 1-entry whose column is labeled as a row. */
    size_t pivotColumn = SIZE_MAX;
    for (size_t column = 0; column < numColumns; ++column)
    {
      if (dense[pivotRow][column] && CMRelementIsRow(columnElements[column]))
      {
        pivotColumn = column;
        break;
      }
    }

    assert(pivotColumn < SIZE_MAX);

    /* Pivot. */
    for (size_t row = 0; row < numRows; ++row)
    {
      if (!dense[row][pivotColumn] || row == pivotRow)
        continue;

      for (size_t column = 0; column < numColumns; ++column)
      {
        if (!dense[pivotRow][column] || column == pivotColumn)
          continue;

        dense[row][column] = dense[row][column] ? 0 : 1;
      }
    }
    CMR_ELEMENT swapElement = rowElements[pivotRow];
    rowElements[pivotRow] = columnElements[pivotColumn];
    columnElements[pivotColumn] = swapElement;
  }

  /* Create mapping from matrix to dense. */
  size_t* rowToDenseRow = (size_t*) malloc( numRows * sizeof(size_t) );
  for (size_t denseRow = 0; denseRow < numRows; ++denseRow)
  {
    assert(CMRelementIsRow(rowElements[denseRow]));
    size_t row = CMRelementToRowIndex(rowElements[denseRow]);
    rowToDenseRow[row] = denseRow;
  }
  size_t* columnToDenseColumn = (size_t*) malloc( numColumns * sizeof(size_t) );
  for (size_t denseColumn = 0; denseColumn < numColumns; ++denseColumn)
  {
    assert(CMRelementIsColumn(columnElements[denseColumn]));
    size_t column = CMRelementToColumnIndex(columnElements[denseColumn]);
    columnToDenseColumn[column] = denseColumn;
  }

  size_t countNonzeros = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    for (size_t column = 0; column < numColumns; ++column)
    {
      if (dense[row][column])
        ++countNonzeros;
    }
  }

  char* message = NULL;
  for (size_t row = 0; row < numRows; ++row)
  {
    size_t first = dec->matrix->rowSlice[row];
    size_t beyond = dec->matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = dec->matrix->entryColumns[e];
      if (!dense[rowToDenseRow[row]][columnToDenseColumn[column]])
      {
        message = CMRconsistencyMessage(
          "nested minor matrix has zero entry at %ld,%ld after pivoting, while the matrix entry is nonzero.",
          row, column);
        break;
      }
    }
    if (message)
      break;
  }

  if (!message && countNonzeros != dec->matrix->numNonzeros)
  {
    message = CMRconsistencyMessage(
      "nested minor matrix has %ld nonzeros after pivoting back, while matrix has %ld nonzeros.",
      countNonzeros, dec->matrix->numNonzeros);
  }

  free(columnToDenseColumn);
  free(rowToDenseRow);
  free(columnElements);
  free(rowElements);
  for (size_t row = 0; row < numRows; ++row)
    free(dense[row]);
  free(dense);

  return message;
}

char* CMRdecConsistency(CMR_DEC* dec, bool recurse)
{
  if (!dec)
    return NULL;

  char* message = consistencyNested(dec);
  if (message)
    return message;

  if (recurse)
  {
    for (size_t c = 0; c < dec->numChildren; ++c)
    {
      char* message = CMRdecConsistency(dec->children[c], true);
      if (message)
      {
        size_t length = strlen(message);
        char* newMessage = (char*) malloc( (length + 16) * sizeof(char));
        sprintf(newMessage, "child %ld: %s", c, newMessage);
        free(message);
        return newMessage;
      }
    }
  }

  return NULL;
}

CMR_GRAPH* CMRdecGraph(CMR_DEC* dec)
{
  assert(dec);

  return dec->graph;
}

CMR_GRAPH_EDGE* CMRdecGraphForest(CMR_DEC* dec)
{
  assert(dec);

  return dec->graphForest;
}

size_t CMRdecGraphSizeForest(CMR_DEC* dec)
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

CMR_GRAPH_EDGE* CMRdecGraphCoforest(CMR_DEC* dec)
{
  assert(dec);

  return dec->graphCoforest;
}

size_t CMRdecGraphSizeCoforest(CMR_DEC* dec)
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

bool* CMRdecGraphArcsReversed(CMR_DEC* dec)
{
  assert(dec);

  return dec->graphArcsReversed;
}

CMR_GRAPH* CMRdecCograph(CMR_DEC* dec)
{
  assert(dec);

  return dec->cograph;
}

CMR_GRAPH_EDGE* CMRdecCographForest(CMR_DEC* dec)
{
  assert(dec);

  return dec->cographForest;
}

CMR_GRAPH_EDGE* CMRdecCographCoforest(CMR_DEC* dec)
{
  assert(dec);

  return dec->cographCoforest;
}

bool* CMRdecCographArcsReversed(CMR_DEC* dec)
{
  assert(dec);

  return dec->cographArcsReversed;
}

size_t CMRdecCographSizeForest(CMR_DEC* dec)
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

size_t CMRdecCographSizeCoforest(CMR_DEC* dec)
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

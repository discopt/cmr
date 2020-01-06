#include <tu/matroid.h>

#include <assert.h>
#include <stdlib.h>

void TUfreeDec(TU* tu, TU_DEC** dec)
{
  assert(dec);
  assert(*dec);
  if ((*dec)->numChildren > 0)
  {
    for (int c = 0; c < (*dec)->numChildren; ++c)
      TUfreeDec(tu, &(*dec)->children[c]);
    TUfreeBlockArray(tu, &(*dec)->children);
  }
  if ((*dec)->matrix)
    TUfreeCharMatrix(tu, &(*dec)->matrix);
  if ((*dec)->transpose)
    TUfreeCharMatrix(tu, &(*dec)->transpose);
  if ((*dec)->graph)
    TUfreeGraph(tu, &(*dec)->graph);
  if ((*dec)->cograph)
    TUfreeGraph(tu, &(*dec)->cograph);
  TUfreeBlock(tu, dec);
}

bool TUisDecLeaf(TU_DEC* dec)
{
  assert(dec);
  return dec->numChildren == 0;
}

bool TUisDecRegular(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_REGULAR;
}

bool TUisDecGraphic(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_GRAPHIC;
}

bool TUisDecCographic(TU_DEC* dec)
{
  assert(dec);
  return dec->flags & TU_DEC_COGRAPHIC;
}

char TUisDecSum(TU_DEC* dec)
{
  assert(dec);
  char result = dec->flags & TU_DEC_TYPE_MASK;
  return result <= 3 ? result : 0;
}

int TUgetDecNumRows(TU_DEC* dec)
{
  assert(dec);
  assert(dec->matrix);
  return dec->matrix->numRows;
}

int TUgetDecNumColumns(TU_DEC* dec)
{
  assert(dec);
  assert(dec->matrix);
  return dec->matrix->numColumns;
}

/**
 * \brief Returns the rank of the lower-left submatrix of this node of the tree.
 *
 * Only valid if \c type is \ref TU_DEC_TWO_SUM or \ref TU_DEC_THREE_SUM.
 */
TU_EXPORT
int TUgetDecRankLowerLeft(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the rank of the top-right submatrix of this node of the tree.
 *
 * Only valid if \c type is \ref TU_DEC_TWO_SUM or \ref TU_DEC_THREE_SUM.
 */
TU_EXPORT
int TUgetDecRankTopRight(
  TU_DEC* dec /**< Decomposition tree */
);

#include <tu/matroid.h>

#include <assert.h>
#include <stdlib.h>

void TUcreateSubmatrix(
  TU_SUBMATRIX** submatrix,
  int numRows,
  int numColumns
  )
{
  assert(submatrix != NULL);
  *submatrix = (TU_SUBMATRIX*) malloc(sizeof(TU_SUBMATRIX));
  (*submatrix)->numRows = numRows;
  (*submatrix)->rows = (int*) malloc(numRows * sizeof(int));
  (*submatrix)->numColumns = numColumns;
  (*submatrix)->columns = (int*) malloc(numColumns * sizeof(int));
}

void TUcreateSubmatrix1x1(
  TU_SUBMATRIX** submatrix, /**< Pointer to submatrix */
  int row, /**< Row of entry */
  int column /**< Column of entry */
  )
{
  TUcreateSubmatrix(submatrix, 1, 1);
  (*submatrix)->rows[0] = row;
  (*submatrix)->columns[0] = column;
}


/**
 * \brief Frees a submatrix.
 */
TU_EXPORT
void TUfreeSubmatrix(
  TU_SUBMATRIX** submatrix /**< Pointer to submatrix */
);

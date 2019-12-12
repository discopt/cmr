#ifndef TU_MATROID_H
#define TU_MATROID_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/env.h>

/**
 * \brief Row and column indices for a submatrix
 * 
 * Does not contain information about the matrix it refers to.
 */
struct TU_Submatrix
{
  /**
   * \brief Number of rows
   */
  int numRows;
  /**
   * \brief Array with row indices
   */
  int* rows;
  /**
   * \brief Number of columns
   */
  int numColumns;
  /**
   * \brief Array with column indices
   */
  int* columns;
};
typedef struct TU_Submatrix TU_SUBMATRIX;

/**
 * \brief Creates a submatrix of given size.
 * 
 * Only allocates the memory. Use \ref TUgetSubmatrixRows and \ref TUgetSubmatrixColumns to modify
 * the row and column indices, respectively.
 */
TU_EXPORT
void TUcreateSubmatrix(
  TU_SUBMATRIX** submatrix, /** Pointer to submatrix */
  int numRows, /**< Number of rows */
  int numColumns /**< Number of columns */
);

/**
 * \brief Creates a 1x1 submatrix.
 */

TU_EXPORT
void TUcreateSubmatrix1x1(
  TU_SUBMATRIX** submatrix, /**< Pointer to submatrix */
  int row, /**< Row of entry */
  int column /**< Column of entry */
);

/**
 * \brief Frees a submatrix.
 */
TU_EXPORT
void TUfreeSubmatrix(
  TU_SUBMATRIX** submatrix /**< Pointer to submatrix */
);
  
enum TU_Dec_Type
{
  TU_DEC_IRREGULAR = 0,
  TU_DEC_GRAPHIC = 1,
  TU_DEC_COGRAPHIC = 2,
  TU_DEC_R10 = 3,
  TU_DEC_ONE_SUM = 4,
  TU_DEC_TWO_SUM = 5,
  TU_DEC_THREE_SUM = 6
};
typedef enum TU_Dec_Type TU_DEC_TYPE;
  
typedef struct TU_Dec TU_DEC;

/**
 * \brief Frees a decomposition tree.
 */
TU_EXPORT
void TUfreeDec(
  TU_DEC** dec /**< Pointer to decomposition */
);

/**
 * \brief Returns the number of matrix rows associated to this decomposition tree.
 */
TU_EXPORT
int TUgetDecNumRows(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the number of matrix columns associated to this decomposition tree.
 */
TU_EXPORT
int TUgetDecNumColumns(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the matrix rows associated to this decomposition tree.
 * 
 * The row indices are returned in the array \p rows which must have memory for at least \p size
 * indices. The method stores at most \p size indices even if there are more. The true number of
 * indices is returned, but can also be queried via \ref TUgetDecNumRows.
 */
TU_EXPORT
int TUgetDecRows(
  TU_DEC* dec, /**< Decomposition tree */
  int* rows, /**< Array for indices */
  int size /**< Length of \p rows array. */
);

/**
 * \brief Returns the matrix columns associated to this decomposition tree.
 * 
 * The column indices are returned in the array \p columns which must have memory for at least
 * \p size indices. The method stores at most \p size indices even if there are more. The true
 * number of indices is returned, but can also be queried via \ref TUgetDecNumColumns.
 */
TU_EXPORT
int TUgetDecColumns(
  TU_DEC* dec, /**< Decomposition tree */
  int* columns, /**< Array for column indices */
  int size /**< Length of \p rows array. */
);


#ifdef __cplusplus
}
#endif

#endif /* TU_MATROID_H */

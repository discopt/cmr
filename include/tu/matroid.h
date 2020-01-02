#ifndef TU_MATROID_H
#define TU_MATROID_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/matrix.h>
#include <tu/graph.h>

typedef enum
{
  TU_DEC_IRREGULAR = 0,
  TU_DEC_GRAPHIC = 1,
  TU_DEC_COGRAPHIC = 2,
  TU_DEC_R10 = 3,
  TU_DEC_ONE_SUM = 4,
  TU_DEC_TWO_SUM = 5,
  TU_DEC_THREE_SUM = 6
} TU_DEC_TYPE;

struct TU_DEC_;

typedef struct TU_DEC_ TU_DEC;

/**
 * \brief Frees a decomposition tree.
 */
TU_EXPORT
void TUfreeDec(
  TU_DEC** dec /**< Pointer to decomposition */
);

/**
 * \brief Returns the type of the root of the decomposition tree.
 */
TU_EXPORT
TU_DEC_TYPE TUgetDecType(
  TU_DEC* dec /**< Decomposition tree */
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
 * \brief Returns the sparse matrix associated to this decomposition tree.
 */
TU_EXPORT
int TUgetDecMatrix(
  TU_DEC* dec, /**< Decomposition tree */
  TU_MATRIX_CHAR* matrix /**< Matrix */
);

/**
 * \brief Returns the transpose of the matrix associated to this decomposition tree.
 */
TU_EXPORT
int TUgetDecTranspose(
  TU_DEC* dec, /**< Decomposition tree */
  TU_MATRIX_CHAR* transpose /**< Transpose of matrix */
);

/**
 * \brief Returns the number of child nodes of the root of this decomposition tree.
 */
TU_EXPORT
int TUgetDecNumChildren(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns a child of the root of this decomposition tree.
 * 
 * \p child must be in [0, \ref TUgetDecNumChildren).
 */
TU_EXPORT
int TUgetDecChild(
  TU_DEC* dec, /**< Decomposition tree */
  int child /**< Index of child. */
);

/**
 * \brief Returns the row labels of the matrix of this decomposition tree.
 */
TU_EXPORT
int* TUgetDecRowLabels(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the column labels of the matrix of this decomposition tree.
 */
TU_EXPORT
int* TUgetDecColumnLabels(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the rank of the lower-left submatrix of the root of this decomposition tree.
 */
TU_EXPORT
int TUgetDecRankLowerLeft(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the rank of the top-right submatrix of the root of this decomposition tree.
 */
TU_EXPORT
int TUgetDecRankTopRight(
  TU_DEC* dec /**< Decomposition tree */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_MATROID_H */

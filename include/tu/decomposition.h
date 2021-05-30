#ifndef TU_DECOMPOSITION_H
#define TU_DECOMPOSITION_H

#include <tu/env.h>
#include <tu/element.h>
#include <tu/matrix.h>
#include <tu/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup Matroid Representation matrices of matroids
 *
 * @{
 */

struct _TU_DEC;

typedef struct _TU_DEC TU_DEC;

typedef enum
{
  TU_DEC_IRREGULAR = -1,          /**< Node represents 3-connected irregular minor. */
  TU_DEC_UNKNOWN = 0,             /**< Type of node is not yet determined. */
  TU_DEC_ONE_SUM = 1,             /**< Node represents a 1-sum of minors; arbitrary number of child nodes. */
  TU_DEC_TWO_SUM = 2,             /**< Node represents a 2-sum of minors; two child nodes. */
  TU_DEC_THREE_SUM = 3,           /**< Node represents a 3-sum of minors; two child nodes. */
  TU_DEC_SUBMATRIX = 4,           /**< Node represents the restriction to a minor; one child node. */
  TU_DEC_PIVOTS = 5,              /**< Node represents a sequence of pivots; one child node. */
  TU_DEC_GRAPHIC = 6,             /**< Node represents a graphic leaf minor; no child nodes. */
  TU_DEC_COGRAPHIC = 7,           /**< Node represents a cographic leaf minor; no child nodes. */
  TU_DEC_PLANAR = 8,              /**< Node represents a planar (graphic and cographic) leaf minor; no child nodes. */

  TU_DEC_SPECIAL_R10 = 16,        /**< Node represents a minor isomorphic to \f$ R_{10} \f$. */
  TU_DEC_SPECIAL_FANO = 17,       /**< Node represents a minor isomorphic to \f$ F_7 \f$. */
  TU_DEC_SPECIAL_FANO_DUAL = 18,  /**< Node represents a minor isomorphic to \f$ F_7^\star \f$. */
  TU_DEC_SPECIAL_K_5 = 19,        /**< Node represents a minor isomorphic to \f$ M(K_5) \f$. */
  TU_DEC_SPECIAL_K_5_DUAL = 20,   /**< Node represents a minor isomorphic to \f$ M(K_5)^\star \f$. */
  TU_DEC_SPECIAL_K_3_3 = 21,      /**< Node represents a minor isomorphic to \f$ M(K_{3,3}) \f$. */
  TU_DEC_SPECIAL_K_3_3_DUAL = 22  /**< Node represents a minor isomorphic to \f$ M(K_{3,3})^\star \f$. */
} TU_DEC_TYPE;

typedef enum
{
  TU_DEC_MASK_REPRESENTATION = 3,       /**< Bit mask for a specific representation matrix of a minor. */
  TU_DEC_IS_GRAPHIC = 4,                /**< Minor is graphic. */
  TU_DEC_IS_COGRAPHIC = 8,              /**< Minor is cographic. */
  TU_DEC_IS_REGULAR = 16,               /**< Minor is regular. */
  TU_DEC_HAS_LOWER_LEFT_NONZEROS = 32,  /**< The 2- or 3-sum has nonzeros in lower-left. */
  TU_DEC_HAS_UPPER_RIGHT_NONZEROS = 64  /**< The 2- or 3-sum has nonzeros in upper-right. */
} TU_DEC_FLAGS;

TU_EXPORT
TU_ERROR TUdecFree(
  TU* tu,       /**< \ref TU environment. */
  TU_DEC** pdec /**< Pointer to the decomposition node. */
);

/**
 * \brief Returns k if \p dec is a k-sum node and 0 otherwise.
 */

TU_EXPORT
int TUdecIsSum(
  TU_DEC* dec,              /**< decomposition. */
  bool* plowerLeftNonzeros, /**< Pointer for storing the lower-left submatrix of a 2- or 3-sum contains a nonzero (may be \c NULL). */
  bool* pupperRightNonzeros /**< Pointer for storing the lower-left submatrix of a 2- or 3-sum contains a nonzero (may be \c NULL). */
);

/**
 * \brief Returns \c true if and only if \p dec is a pivot sequence node.
 */

TU_EXPORT
bool TUdecIsPivotSequence(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is a submatrix node.
 */

TU_EXPORT
bool TUdecIsSubmatrix(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the corresponding flag if \p dec is a special matrix node and 0 otherwise.
 */

TU_EXPORT
TU_DEC_TYPE TUdecIsSpecialLeaf(
  TU_DEC* dec,                /**< decomposition. */
  int* prepresentationMatrix  /**< Pointer for storing the id of the actual representation matrix (may be \c NULL). */
);

/**
 * \brief Returns \c true if and only if \p dec is a graphic leaf node.
 */

TU_EXPORT
bool TUdecIsGraphicLeaf(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is a cographic leaf node.
 */

TU_EXPORT
bool TUdecIsCographicLeaf(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is graphic.
 */

TU_EXPORT
bool TUdecIsGraphic(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is cographic.
 */

TU_EXPORT
bool TUdecIsCographic(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is regular.
 */

TU_EXPORT
bool TUdecIsRegular(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the number of rows.
 */

TU_EXPORT
bool TUdecNumRows(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the row elements.
 */

TU_EXPORT
TU_ELEMENT* TUdecRowElements(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the mapping of rows to rows of parent.
 */

TU_EXPORT
size_t* TUdecRowsParent(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the number of columns.
 */

TU_EXPORT
bool TUdecNumColumns(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the column elements.
 */

TU_EXPORT
TU_ELEMENT* TUdecColumnElements(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the mapping of columns to columns of parent.
 */

TU_EXPORT
size_t* TUdecColumnsParent(
  TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Prints the decomposition \p dec to \p stream.
 */

TU_EXPORT
TU_ERROR TUdecPrint(
  FILE* stream, /**< Stream. */
  TU_DEC* dec,  /**< Decomposition node. */
  size_t indent /**< Indentation of this node. */
);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* TU_DECOMPOSITION_H */

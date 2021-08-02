#ifndef CMR_TU_DECOMPOSITION_H
#define CMR_TU_DECOMPOSITION_H

#include <cmr/env.h>
#include <cmr/element.h>
#include <cmr/matrix.h>
#include <cmr/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup Matroid Representation matrices of matroids
 *
 * @{
 */

struct _CMR_TU_DEC;

typedef struct _CMR_TU_DEC CMR_TU_DEC;

typedef enum
{
  CMR_TU_DEC_IRREGULAR = -1,          /**< Node represents 3-connected irregular minor. */
  CMR_TU_DEC_UNKNOWN = 0,             /**< Type of node is not yet determined. */
  CMR_TU_DEC_ONE_SUM = 1,             /**< Node represents a 1-sum of minors; arbitrary number of child nodes. */
  CMR_TU_DEC_TWO_SUM = 2,             /**< Node represents a 2-sum of minors; two child nodes. */
  CMR_TU_DEC_THREE_SUM = 3,           /**< Node represents a 3-sum of minors; two child nodes. */
  CMR_TU_DEC_SUBMATRIX = 4,           /**< Node represents the restriction to a minor; one child node. */
  CMR_TU_DEC_PIVOTS = 5,              /**< Node represents a sequence of pivots; one child node. */
  CMR_TU_DEC_GRAPHIC = 6,             /**< Node represents a graphic leaf minor; no child nodes. */
  CMR_TU_DEC_COGRAPHIC = 7,           /**< Node represents a cographic leaf minor; no child nodes. */
  CMR_TU_DEC_PLANAR = 8,              /**< Node represents a planar (graphic and cographic) leaf minor; no child nodes. */

  CMR_TU_DEC_SPECIAL_R10 = 16,        /**< Node represents a minor isomorphic to \f$ R_{10} \f$. */
  CMR_TU_DEC_SPECIAL_FANO = 17,       /**< Node represents a minor isomorphic to \f$ F_7 \f$. */
  CMR_TU_DEC_SPECIAL_FANO_DUAL = 18,  /**< Node represents a minor isomorphic to \f$ F_7^\star \f$. */
  CMR_TU_DEC_SPECIAL_K_5 = 19,        /**< Node represents a minor isomorphic to \f$ M(K_5) \f$. */
  CMR_TU_DEC_SPECIAL_K_5_DUAL = 20,   /**< Node represents a minor isomorphic to \f$ M(K_5)^\star \f$. */
  CMR_TU_DEC_SPECIAL_K_3_3 = 21,      /**< Node represents a minor isomorphic to \f$ M(K_{3,3}) \f$. */
  CMR_TU_DEC_SPECIAL_K_3_3_DUAL = 22  /**< Node represents a minor isomorphic to \f$ M(K_{3,3})^\star \f$. */
} CMR_TU_DEC_TYPE;

typedef enum
{
  CMR_TU_DEC_MASK_REPRESENTATION = 3,       /**< Bit mask for a specific representation matrix of a minor. */
  CMR_TU_DEC_IS_GRAPHIC = 4,                /**< Minor is graphic. */
  CMR_TU_DEC_IS_COGRAPHIC = 8,              /**< Minor is cographic. */
  CMR_TU_DEC_IS_REGULAR = 16,               /**< Minor is regular. */
  CMR_TU_DEC_HAS_LOWER_LEFT_NONZEROS = 32,  /**< The 2- or 3-sum has nonzeros in lower-left. */
  CMR_TU_DEC_HAS_UPPER_RIGHT_NONZEROS = 64  /**< The 2- or 3-sum has nonzeros in upper-right. */
} CMR_TU_DEC_FLAGS;

CMR_EXPORT
CMR_ERROR CMRtudecFree(
  CMR* cmr,       /**< \ref CMR environment. */
  CMR_TU_DEC** pdec /**< Pointer to the decomposition node. */
);

/**
 * \brief Returns k if \p dec is a k-sum node and 0 otherwise.
 */

CMR_EXPORT
int CMRtudecIsSum(
  CMR_TU_DEC* dec,              /**< decomposition. */
  bool* plowerLeftNonzeros, /**< Pointer for storing the lower-left submatrix of a 2- or 3-sum contains a nonzero (may be \c NULL). */
  bool* pupperRightNonzeros /**< Pointer for storing the lower-left submatrix of a 2- or 3-sum contains a nonzero (may be \c NULL). */
);

/**
 * \brief Returns \c true if and only if \p dec is a pivot sequence node.
 */

CMR_EXPORT
bool CMRtudecIsPivotSequence(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is a submatrix node.
 */

CMR_EXPORT
bool CMRtudecIsSubmatrix(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the corresponding flag if \p dec is a special matrix node and 0 otherwise.
 */

CMR_EXPORT
CMR_TU_DEC_TYPE CMRtudecIsSpecialLeaf(
  CMR_TU_DEC* dec,                /**< decomposition. */
  int* prepresentationMatrix  /**< Pointer for storing the id of the actual representation matrix (may be \c NULL). */
);

/**
 * \brief Returns \c true if and only if \p dec is a graphic leaf node.
 */

CMR_EXPORT
bool CMRtudecIsGraphicLeaf(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is a cographic leaf node.
 */

CMR_EXPORT
bool CMRtudecIsCographicLeaf(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is graphic.
 */

CMR_EXPORT
bool CMRtudecIsGraphic(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is cographic.
 */

CMR_EXPORT
bool CMRtudecIsCographic(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is regular.
 */

CMR_EXPORT
bool CMRtudecIsRegular(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the number of rows.
 */

CMR_EXPORT
bool CMRtudecNumRows(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the row elements.
 */

CMR_EXPORT
CMR_ELEMENT* CMRtudecRowElements(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the mapping of rows to rows of parent.
 */

CMR_EXPORT
size_t* CMRtudecRowsParent(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the number of columns.
 */

CMR_EXPORT
bool CMRtudecNumColumns(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the column elements.
 */

CMR_EXPORT
CMR_ELEMENT* CMRtudecColumnElements(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Returns the mapping of columns to columns of parent.
 */

CMR_EXPORT
size_t* CMRtudecColumnsParent(
  CMR_TU_DEC* dec /**< decomposition. */
);

/**
 * \brief Prints the decomposition \p dec to \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRtudecPrint(
  FILE* stream, /**< Stream. */
  CMR_TU_DEC* dec,  /**< Decomposition node. */
  size_t indent /**< Indentation of this node. */
);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_TU_DECOMPOSITION_H */

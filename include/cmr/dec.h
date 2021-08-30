#ifndef CMR_DEC_H
#define CMR_DEC_H

/**
 * \file dec.h
 *
 * \author Matthias Walter
 *
 * \brief Data structures for a decomposition tree for a matrix.
 */

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

struct _CMR_DEC;

typedef struct _CMR_DEC CMR_DEC;

typedef enum
{
  CMR_DEC_IRREGULAR = -1,          /**< Node represents 3-connected irregular minor. */
  CMR_DEC_UNKNOWN = 0,             /**< Type of node is not yet determined. */
  CMR_DEC_ONE_SUM = 1,             /**< Node represents a 1-sum of minors; arbitrary number of child nodes. */
  CMR_DEC_TWO_SUM = 2,             /**< Node represents a 2-sum of minors; two child nodes. */
  CMR_DEC_THREE_SUM = 3,           /**< Node represents a 3-sum of minors; two child nodes. */
  CMR_DEC_GRAPHIC = 4,             /**< Node represents a graphic leaf minor; no child nodes. */
  CMR_DEC_COGRAPHIC = 5,           /**< Node represents a cographic leaf minor; no child nodes. */
  CMR_DEC_PLANAR = 6,              /**< Node represents a planar (graphic and cographic) leaf minor; no child nodes. */

  CMR_DEC_SPECIAL_R10 = 16,        /**< Node represents a minor isomorphic to \f$ R_{10} \f$. */
  CMR_DEC_SPECIAL_FANO = 17,       /**< Node represents a minor isomorphic to \f$ F_7 \f$. */
  CMR_DEC_SPECIAL_FANO_DUAL = 18,  /**< Node represents a minor isomorphic to \f$ F_7^\star \f$. */
  CMR_DEC_SPECIAL_K_5 = 19,        /**< Node represents a minor isomorphic to \f$ M(K_5) \f$. */
  CMR_DEC_SPECIAL_K_5_DUAL = 20,   /**< Node represents a minor isomorphic to \f$ M(K_5)^\star \f$. */
  CMR_DEC_SPECIAL_K_3_3 = 21,      /**< Node represents a minor isomorphic to \f$ M(K_{3,3}) \f$. */
  CMR_DEC_SPECIAL_K_3_3_DUAL = 22  /**< Node represents a minor isomorphic to \f$ M(K_{3,3})^\star \f$. */
} CMR_DEC_TYPE;

typedef enum
{
  CMR_DEC_MASK_REPRESENTATION = 3,       /**< Bit mask for a specific representation matrix of a minor. */
  CMR_DEC_IS_GRAPHIC = 4,                /**< Minor is graphic. */
  CMR_DEC_IS_COGRAPHIC = 8,              /**< Minor is cographic. */
  CMR_DEC_IS_REGULAR = 16,               /**< Minor is regular. */
  CMR_DEC_HAS_LOWER_LEFT_NONZEROS = 32,  /**< The 2- or 3-sum has nonzeros in lower-left. */
  CMR_DEC_HAS_UPPER_RIGHT_NONZEROS = 64  /**< The 2- or 3-sum has nonzeros in upper-right. */
} CMR_DEC_FLAGS;

CMR_EXPORT
CMR_ERROR CMRdecFree(
  CMR* cmr,       /**< \ref CMR environment. */
  CMR_DEC** pdec  /**< Pointer to the decomposition node. */
);

/**
 * \brief Returns k if \p dec is a k-sum node and 0 otherwise.
 */

CMR_EXPORT
int CMRdecIsSum(
  CMR_DEC* dec,             /**< Decomposition node. */
  bool* plowerLeftNonzeros, /**< Pointer for storing the lower-left submatrix of a 2- or 3-sum contains a nonzero (may be \c NULL). */
  bool* pupperRightNonzeros /**< Pointer for storing the lower-left submatrix of a 2- or 3-sum contains a nonzero (may be \c NULL). */
);

/**
 * \brief Returns the corresponding flag if \p dec is a special matrix node and 0 otherwise.
 */

CMR_EXPORT
CMR_DEC_TYPE CMRdecIsSpecialLeaf(
  CMR_DEC* dec,               /**< Decomposition. */
  int* prepresentationMatrix  /**< Pointer for storing the id of the actual representation matrix (may be \c NULL). */
);

/**
 * \brief Returns \c true if and only if \p dec is a graphic leaf node.
 */

CMR_EXPORT
bool CMRdecIsGraphicLeaf(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is a cographic leaf node.
 */

CMR_EXPORT
bool CMRdecIsCographicLeaf(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is graphic.
 */

CMR_EXPORT
bool CMRdecIsGraphic(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is cographic.
 */

CMR_EXPORT
bool CMRdecIsCographic(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns \c true if and only if \p dec is regular.
 */

CMR_EXPORT
bool CMRdecIsRegular(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns the number of rows.
 */

CMR_EXPORT
bool CMRdecNumRows(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns the row elements.
 */

CMR_EXPORT
CMR_ELEMENT* CMRdecRowElements(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns the mapping of rows to rows of parent.
 */

CMR_EXPORT
size_t* CMRdecRowsParent(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns the number of columns.
 */

CMR_EXPORT
bool CMRdecNumColumns(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns the column elements.
 */

CMR_EXPORT
CMR_ELEMENT* CMRdecColumnElements(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Returns the mapping of columns to columns of parent.
 */

CMR_EXPORT
size_t* CMRdecColumnsParent(
  CMR_DEC* dec  /**< Decomposition. */
);

/**
 * \brief Prints the decomposition \p dec to \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRdecPrint(
  FILE* stream, /**< Stream. */
  CMR_DEC* dec, /**< Decomposition. */
  size_t indent /**< Indentation of this node. */
);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_DEC_H */

#ifndef CMR_REGULAR_H
#define CMR_REGULAR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/graph.h>

typedef enum
{
  CMR_TU_DEC_ONE_SUM = 1,       /**< 1-sum decomposition with multiple children. */
  CMR_TU_DEC_TWO_SUM = 2,       /**< 2-sum decomposition with 2 children. */
  CMR_TU_DEC_THREE_SUM = 3,     /**< 3-sum decomposition with 2 children. */
  CMR_TU_DEC_R10 = 4,           /**< equal to R10. */
  CMR_TU_DEC_MULTI_SUM = 5,     /**< Sequence of simple 2- or 3-sums. */
  CMR_TU_DEC_TYPE_MASK = 7,
  CMR_TU_DEC_PROCESSED = 8,     /**< Node was further decomposed or identified as a leaf. */
  CMR_TU_DEC_GRAPHIC = 16,      /**< Node is graphic; also set for sums of graphic matroids. */
  CMR_TU_DEC_COGRAPHIC = 32,    /**< Node is cographic; also set for sums of cographic matroids. */
  CMR_TU_DEC_REGULAR = 64,      /**< Node is regular; also set for sums of regular matroids. */
  CMR_TU_DEC_RANK_LOWER_LEFT = 128, /**< The 2- or 3-sum has nonzeros in lower-left. */
  CMR_TU_DEC_RANK_UPPER_RIGHT = 256 /**< The 2- or 3-sum has nonzeros in upper-right. */
} CMR_TU_DEC_FLAGS;

typedef struct _CMR_TU_DEC
{
  CMR_CHRMAT* matrix;     /**< Binary matrix representing this tree node's matroid. */
  CMR_CHRMAT* transpose;  /**< Transpose of matrix representing this tree node's matroid. */
  int* rowLabels;             /**< Array with row labels. */
  int* columnLabels;          /**< Array with column labels. */
  int* parentRows;            /**< Mapping to rows of parent node or \c NULL. */
  int* parentColumns;         /**< Mapping to columns of parent node or \c NULL. */

  CMR_TU_DEC_FLAGS flags;         /**< Flags for this tree node. */
  CMR_GRAPH *graph;        /**< Graph corresponding to this tree node, or \c NULL.  */
  CMR_GRAPH *cograph;      /**< Cograph corresponding to this tree node, or \c NULL.  */
  int numChildren;            /**< Number of children of this tree node. */
  struct _CMR_TU_DEC** children;  /**< Array with pointers to children of this tree node. */
} CMR_TU_DEC;

/**
 * \brief Frees a decomposition tree (including the data of all its nodes.
 */
CMR_EXPORT
void CMRtudecFree(
  CMR* cmr,       /**< \ref CMR environment. */
  CMR_TU_DEC** pdec /**< Pointer to decomposition tree. */
);

/**
 * \brief Returns \c true iff this node is a leaf of the tree.
 */
CMR_EXPORT
bool CMRtudecIsLeaf(
  CMR_TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns \c true if and only if tree node corresponds to a regular matroid.
 */
CMR_EXPORT
bool CMRtudecIsRegular(
  CMR_TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns \c true if and only if tree node corresponds to a graphic matroid.
 */
CMR_EXPORT
bool CMRtudecIsGraphic(
  CMR_TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns \c true if and only if tree node corresponds to a cographic matroid.
 */
CMR_EXPORT
bool CMRtudecIsCographic(
  CMR_TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns \c k if tree node corresponds to a k-decomposition, and 0 otherwise.
 */
CMR_EXPORT
char CMRtudecIsSum(
  CMR_TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the number of matrix rows associated to this decomposition tree.
 */
CMR_EXPORT
int CMRtudecNumRows(
  CMR_TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the number of matrix columns associated to this decomposition tree.
 */
CMR_EXPORT
int CMRtudecNumColumns(
  CMR_TU_DEC* dec /**< Decomposition tree */
);


/**
 * \brief Returns the rank of the lower-left submatrix of this node of the tree.
 *
 * Only valid if either flag \ref CMR_DEC_TWO_SUM or \ref CMR_DEC_THREE_SUM is set.
 */
CMR_EXPORT
int CMRgetDecRankLowerLeft(
  CMR_TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the rank of the upper-right submatrix of this node of the tree.
 *
 * Only valid if either flag \ref CMR_DEC_TWO_SUM or \ref CMR_DEC_THREE_SUM is set.
 */
CMR_EXPORT
int CMRgetDecRankUpperRight(
  CMR_TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Tests binary linear matroid for regularity.
 *
 * If \p decomposition is not \c NULL, \c *decomposition will be a decomposition tree.
 *
 */

CMR_EXPORT
bool CMRregularTest(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,  /**< Char matrix. */
  int* rowLabels,     /**< Row labels of matrix (may be \c NULL). */
  int* columnLabels,  /**< Column labels of matrix (may be \c NULL). */
  CMR_TU_DEC** pdec       /**< Pointer for storing the decomposition tree (may be \c NULL). */
);

/**
 * \brief Checks for a 1-sum decomposition of a given ternary or binary linear matroid.
 *
 * Initializes the \c *pdecomposition to a partial decomposition tree. If the matroid is a 1-sum
 * then \c *pdecomposition will be a \ref CMR_DEC_ONE_SUM tree node and its children will be
 * initialized with corresponding sequentially connected submatrices. Otherwise, \c *pdecomposition
 * itself will be initialized with a sequentially connected permutation of the given matroid.
 */

CMR_EXPORT
int CMRregularDecomposeOneSum(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,          /**< Given matrix. */
  int* rowLabels,             /**< Row labels of matrix (may be \c NULL). */
  int* columnLabels,          /**< Column labels of matrix (may be \c NULL). */
  CMR_TU_DEC** pdec,              /**< Pointer for storing the decomposition tree (may be \c NULL). */
  bool constructDecomposition /**< Whether to construct the decomposition. */
);

/**
 * \brief Checks for a simple 2-sum or 3-sum decomposition of a given ternary or binary linear
 *        matroid.
 *
 * Assumes that \p decomposition contains valid \c matrix and \c transpose entries.
 * A simple 2-sum is one where one component consists of only one row and one column. It is found
 * by splitting off unit vectors.
 *
 * \returns Remaining submatrix.
 */

CMR_EXPORT
CMR_TU_DEC* CMRregularDecomposeSimpleSums(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_TU_DEC* dec,                /**< The partial decomposition that is produced. */
  bool unitVectors,           /**< Split off unit vectors via 2-sums. */
  bool paths,                 /**< Split off paths via 3-sums. */
  bool constructDecomposition /**< Whether to construct a proper decomposition. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_REGULAR_H */


#ifndef TU_REGULAR_H
#define TU_REGULAR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/env.h>
#include <tu/matrix.h>
#include <tu/graph.h>

typedef enum
{
  TU_DEC_ONE_SUM = 1,       /**< 1-sum decomposition with multiple children. */
  TU_DEC_TWO_SUM = 2,       /**< 2-sum decomposition with 2 children. */
  TU_DEC_THREE_SUM = 3,     /**< 3-sum decomposition with 2 children. */
  TU_DEC_R10 = 4,           /**< equal to R10. */
  TU_DEC_MULTI_SUM = 5,     /**< Sequence of simple 2- or 3-sums. */
  TU_DEC_TYPE_MASK = 7,
  TU_DEC_PROCESSED = 8,     /**< Node was further decomposed or identified as a leaf. */
  TU_DEC_GRAPHIC = 16,      /**< Node is graphic; also set for sums of graphic matroids. */
  TU_DEC_COGRAPHIC = 32,    /**< Node is cographic; also set for sums of cographic matroids. */
  TU_DEC_REGULAR = 64,      /**< Node is regular; also set for sums of regular matroids. */
  TU_DEC_RANK_LOWER_LEFT = 128, /**< The 2- or 3-sum has nonzeros in lower-left. */
  TU_DEC_RANK_UPPER_RIGHT = 256 /**< The 2- or 3-sum has nonzeros in upper-right. */
} TU_DEC_FLAGS;

typedef struct _TU_DEC
{
  TU_CHRMAT* matrix;     /**< Binary matrix representing this tree node's matroid. */
  TU_CHRMAT* transpose;  /**< Transpose of matrix representing this tree node's matroid. */
  int* rowLabels;             /**< Array with row labels. */
  int* columnLabels;          /**< Array with column labels. */
  int* parentRows;            /**< Mapping to rows of parent node or \c NULL. */
  int* parentColumns;         /**< Mapping to columns of parent node or \c NULL. */

  TU_DEC_FLAGS flags;         /**< Flags for this tree node. */
  TU_GRAPH *graph;        /**< Graph corresponding to this tree node, or \c NULL.  */
  TU_GRAPH *cograph;      /**< Cograph corresponding to this tree node, or \c NULL.  */
  int numChildren;            /**< Number of children of this tree node. */
  struct _TU_DEC** children;  /**< Array with pointers to children of this tree node. */
} TU_DEC;

/**
 * \brief Frees a decomposition tree (including the data of all its nodes.
 */
TU_EXPORT
void TUdecFree(
  TU* tu,       /**< TU environment. */
  TU_DEC** dec  /**< Pointer to decomposition tree. */
);

/**
 * \brief Returns \c true iff this node is a leaf of the tree.
 */
TU_EXPORT
bool TUdecIsLeaf(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns \c true if and only if tree node corresponds to a regular matroid.
 */
TU_EXPORT
bool TUdecIsRegular(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns \c true if and only if tree node corresponds to a graphic matroid.
 */
TU_EXPORT
bool TUdecIsGraphic(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns \c true if and only if tree node corresponds to a cographic matroid.
 */
TU_EXPORT
bool TUdecIsCographic(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns \c k if tree node corresponds to a k-decomposition, and 0 otherwise.
 */
TU_EXPORT
char TUdecIsSum(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the number of matrix rows associated to this decomposition tree.
 */
TU_EXPORT
int TUdecNumRows(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the number of matrix columns associated to this decomposition tree.
 */
TU_EXPORT
int TUdecNumColumns(
  TU_DEC* dec /**< Decomposition tree */
);


/**
 * \brief Returns the rank of the lower-left submatrix of this node of the tree.
 *
 * Only valid if either flag \ref TU_DEC_TWO_SUM or \ref TU_DEC_THREE_SUM is set.
 */
TU_EXPORT
int TUgetDecRankLowerLeft(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Returns the rank of the upper-right submatrix of this node of the tree.
 *
 * Only valid if either flag \ref TU_DEC_TWO_SUM or \ref TU_DEC_THREE_SUM is set.
 */
TU_EXPORT
int TUgetDecRankUpperRight(
  TU_DEC* dec /**< Decomposition tree */
);

/**
 * \brief Tests binary linear matroid for regularity.
 *
 * If \p decomposition is not \c NULL, \c *decomposition will be a decomposition tree.
 *
 */

TU_EXPORT
bool TUregularTest(
  TU* tu,                 /**< TU environment */
  TU_CHRMAT* matrix, /**< Char matrix */
  int* rowLabels,         /**< Row labels of matrix; can be \c NULL. */
  int* columnLabels,      /**< Column labels of matrix; can be \c NULL. */
  TU_DEC** decomposition  /**< If not \c NULL, the decomposition tree is stored. */
);

/**
 * \brief Checks for a 1-sum decomposition of a given ternary or binary linear matroid.
 *
 * Initializes the \c *pdecomposition to a partial decomposition tree. If the matroid is a 1-sum
 * then \c *pdecomposition will be a \ref TU_DEC_ONE_SUM tree node and its children will be
 * initialized with corresponding sequentially connected submatrices. Otherwise, \c *pdecomposition
 * itself will be initialized with a sequentially connected permutation of the given matroid.
 */

TU_EXPORT
int TUregularDecomposeOneSum(
  TU* tu,                     /**< TU environment. */
  TU_CHRMAT* matrix,     /**< Given matrix. */
  int* rowLabels,             /**< Row labels of matrix; can be \c NULL. */
  int* columnLabels,          /**< Column labels of matrix; can be \c NULL. */
  TU_DEC** pdecomposition,    /**< Pointer for storing the partial decomposition. */
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

TU_EXPORT
TU_DEC* TUregularDecomposeSimpleSums(
  TU* tu,                     /**< TU environment. */
  TU_DEC* decomposition,      /**< The partial decomposition that is produced. */
  bool unitVectors,           /**< Split off unit vectors via 2-sums. */
  bool paths,                 /**< Split off paths via 3-sums. */
  bool constructDecomposition /**< Whether to construct a proper decomposition. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_REGULAR_H */


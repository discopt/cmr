#ifndef CMR_SEYMOUR_H
#define CMR_SEYMOUR_H

/**
 * \file seymour.h
 *
 * \author Matthias Walter
 *
 * \brief Functionality for Seymour decomposition.
 */

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/matroid.h>
#include <cmr/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup Seymour Seymour decomposition
 *
 * @{
 */

struct _CMR_SEYMOUR_NODE;

typedef struct _CMR_SEYMOUR_NODE CMR_SEYMOUR_NODE;

typedef enum
{
  CMR_SEYMOUR_NODE_TYPE_IRREGULAR = -1,
    /**< Node represents 3-connected irregular matrix. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_UNKNOWN = 0,
    /**< Type of node is not yet determined. */
  CMR_SEYMOUR_NODE_TYPE_ONE_SUM = 1,
    /**< Node represents a \f$ 1 \f$-sum of matrices; has at least 2 child nodes. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_TWO_SUM = 2,
    /**< Node represents a 2-sum of matrices; has two child nodes. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_THREE_SUM = 3,
    /**< Node represents a 3-sum of matrices; has two child nodes. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL = 4,
    /**< Node represents a series-parallel reduction; has one child node. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_PIVOTS = 5,
    /**< Node represents an application of pivots; has one child node. \see \ref seymour_decomposition. */

  CMR_SEYMOUR_NODE_TYPE_GRAPH = 6,
    /**< Node represents a graphic or network leaf. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_COGRAPH = 7,
    /**< Node represents a cographic or conetwork leaf. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_PLANAR = 8,
    /**< Node represents a graphic and cographic (network and conetwork) leaf.  \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_R10 = 9,
    /**< Node represents a representation matrix of \f$ R_{10} \f$.  \see \ref seymour_decomposition. */
} CMR_SEYMOUR_NODE_TYPE;

/**
 * \brief Flags that indicate the type of \f$ 3 \f$-separation.
 *
 * \see The desired types can be set by modifying \ref CMR_REGULAR_PARAMS.threeSumStrategy.
 **/

typedef enum
{
  CMR_SEYMOUR_NODE_THREESUM_FLAG_NO_PIVOTS = 0,
    /**< Indicate to not change the rank distribution; only serves as an option; each constructed node will have either
     **  \ref CMR_SEYMOUR_NODE_THREESUM_FLAG_DISTRIBUTED_RANKS or
     **  \ref CMR_SEYMOUR_NODE_THREESUM_FLAG_CONCENTRATED_RANK set. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_DISTRIBUTED_RANKS = 1,
    /**< The two off-diagonal submatrices both have rank 1.
     **  \see \ref matroid_decomposition. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_CONCENTRATED_RANK = 2,
    /**< The bottom-left submatrix has rank 2 and the top-right submatrix has rank 0.
     **  \see \ref matroid_decomposition. */

  CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_WIDE = 4,
    /**< The first child node is of the form \f$ M_1^{\text{wide}} \f$; valid for distributed ranks.
     **  \see \ref matroid_decomposition. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_TALL = 8,
    /**< The first child node is of the form \f$ M_1^{\text{tall}} \f$; valid for distributed ranks.
     **  \see \ref matroid_decomposition. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_MIXED = 64,
    /**< The first child node is of the form \f$ M_1^{\text{mixed}} \f$; valid for concentrated rank.
     **  \see \ref matroid_decomposition. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_ALLREPR = 128,
    /**< The first child node is of the form \f$ M_1^{\text{all-repr}} \f$; valid for concentrated rank.
     **  \see \ref matroid_decomposition. */

  CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_WIDE = 16,
    /**< The second child node is of the form \f$ M_2^{\text{wide}} \f$; valid for distributed ranks.
     **  \see \ref matroid_decomposition. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_TALL = 32,
    /**< The second child node is of the form \f$ M_2^{\text{tall}} \f$; valid for distributed ranks.
     **  \see \ref matroid_decomposition. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_MIXED = 256,
    /**< The second child node is of the form \f$ M_2^{\text{mixed}} \f$; valid for concentrated rank.
     **  \see \ref matroid_decomposition. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_ALLREPR = 512,
    /**< The second child node is of the form \f$ M_2^{\text{all-repr}} \f$; valid for concentrated rank.
     **  \see \ref matroid_decomposition. */

  CMR_SEYMOUR_NODE_THREESUM_FLAG_SEYMOUR = CMR_SEYMOUR_NODE_THREESUM_FLAG_DISTRIBUTED_RANKS
    | CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_WIDE | CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_WIDE,
    /**< This combination of flags indicates a \f$ 3 \f$-sum as defined by Seymour.
     **  \see \ref matroid_decomposition. */
  CMR_SEYMOUR_NODE_THREESUM_FLAG_TRUEMPER = CMR_SEYMOUR_NODE_THREESUM_FLAG_CONCENTRATED_RANK
    | CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_MIXED | CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_MIXED,
    /**< This combination of flags indicates a \f$ 3 \f$-sum as defined by Truemper.
     **  \see \ref matroid_decomposition. */
} CMR_SEYMOUR_NODE_THREESUM_FLAG;

/**
 * \brief Returns \c true iff the decomposition is over \f$ \mathbb{F}_3 \f$.
 */

CMR_EXPORT
bool CMRseymourIsTernary(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns \c true iff the 3-sum decomposition node has \f$ 3 \f$-separation with two rank-1 matrices.
 */

CMR_EXPORT
bool CMRseymourThreeSumDistributedRanks(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns \c true iff the 3-sum decomposition node has \f$ 3 \f$-separation with one rank-2 matrix.
 */

CMR_EXPORT
bool CMRseymourThreeSumConcentratedRank(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns \c true iff the transposed matrix of the decomposition node \p dec is stored.
 */

CMR_EXPORT
bool CMRseymourHasTranspose(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the matrix of the decomposition node \p dec (or \c NULL if it is not stored).
 */

CMR_EXPORT
CMR_CHRMAT* CMRseymourGetMatrix(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the transposed matrix of the decomposition node \p dec (or \c NULL if it is not stored).
 */

CMR_EXPORT
CMR_CHRMAT* CMRseymourGetTranspose(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of children of the decomposition node \p dec.
 */

CMR_EXPORT
size_t CMRseymourNumChildren(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns a child of the decomposition node \p dec.
 */

CMR_EXPORT
CMR_SEYMOUR_NODE* CMRseymourChild(
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  size_t childIndex         /**< Index of child. */
);

/**
 * \brief Returns the type of a decomposition node \p dec.
 */

CMR_EXPORT
CMR_SEYMOUR_NODE_TYPE CMRseymourType(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of minors of the decomposition node.
 */

CMR_EXPORT
size_t CMRseymourNumMinors(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns a minor of the decomposition node.
 */

CMR_EXPORT
CMR_MINOR* CMRseymourMinor(
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  size_t minorIndex         /**< Index of minor. */
);

/**
 * \brief Indicates graphicness/being network.
 *
 * Returns a positive value if the matrix corresponding to \p dec is graphic/network, zero if it is not known
 * and a negative value otherwise.
 */

CMR_EXPORT
int8_t CMRseymourGraphicness(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Indicates cographicness/being conetwork.
 *
 * Returns a positive value if the matrix corresponding to \p dec is cographic/conetwork, zero if it is not known
 * and a negative value otherwise.
 */

CMR_EXPORT
int8_t CMRseymourCographicness(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Indicates regularity/total unimodularity.
 *
 * Returns a positive value if the matrix corresponding to \p dec is regular/totally unimodular, zero if it is not known
 * and a negative value otherwise.
 */

CMR_EXPORT
int8_t CMRseymourRegularity(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of rows.
 */

CMR_EXPORT
size_t CMRseymourNumRows(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of columns.
 */

CMR_EXPORT
size_t CMRseymourNumColumns(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the mapping of rows of child \p childIndex to this node's elements.
 */

CMR_EXPORT
CMR_ELEMENT* CMRseymourChildRowsToParent(
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  size_t childIndex       /**< Index of child to consider. */
);

/**
 * \brief Returns the mapping of columns of child \p childIndex to this node's elements.
 */

CMR_EXPORT
CMR_ELEMENT* CMRseymourChildColumnsToParent(
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  size_t childIndex       /**< Index of child to consider. */
);

/**
 * \brief Returns the graph (if available).
 */

CMR_EXPORT
CMR_GRAPH* CMRseymourGraph(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the forest of the graph (if available).
 */

CMR_EXPORT
CMR_GRAPH_EDGE* CMRseymourGraphForest(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of edges of the graph's forest (if available).
 */

CMR_EXPORT
size_t CMRseymourGraphSizeForest(
    CMR_SEYMOUR_NODE* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the coforest of the graph (if available).
 */

CMR_EXPORT
CMR_GRAPH_EDGE* CMRseymourGraphCoforest(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of edges of the graph's coforest (if available).
 */

CMR_EXPORT
size_t CMRseymourGraphSizeCoforest(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns an array that indicates for the graph's edges whether they must be reversed (if available).
 */

CMR_EXPORT
bool* CMRseymourGraphArcsReversed(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the cograph (if available).
 */

CMR_EXPORT
CMR_GRAPH* CMRseymourCograph(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of edges of the cograph's forest (if available).
 */

CMR_EXPORT
size_t CMRseymourCographSizeForest(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the forest of the cograph (if available).
 */

CMR_EXPORT
CMR_GRAPH_EDGE* CMRseymourCographForest(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of edges of the cograph's coforest (if available).
 */

CMR_EXPORT
size_t CMRseymourCographSizeCoforest(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the coforest of the cograph (if available).
 */

CMR_EXPORT
CMR_GRAPH_EDGE* CMRseymourCographCoforest(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns an array that indicates for the cograph's edges whether they must be reversed (if available).
 */

CMR_EXPORT
bool* CMRseymourCographArcsReversed(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the number of pivots (if available).
 */

CMR_EXPORT
size_t CMRseymourNumPivots(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the array with the pivot rows (if available).
 */

CMR_EXPORT
size_t* CMRseymourPivotRows(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Returns the array with the pivot columns (if available).
 */

CMR_EXPORT
size_t* CMRseymourPivotColumns(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Prints the decomposition \p dec to \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRseymourPrint(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node,   /**< Seymour decomposition node. */
  FILE* stream,             /**< Stream to write to. */
  bool printChildren,       /**< Whether to recurse. */
  bool printParentElements, /**< Whether to print mapping of rows/columns to parent elements. */
  bool printMatrices,       /**< Whether to print matrices. */
  bool printGraphs,         /**< Whether to print graphs. */
  bool printReductions,     /**< Whether to print series-parallel reductions. */
  bool printPivots          /**< Whether to print pivots. */
);

/**
 * \brief Increases the reference counter by 1.
 */

CMR_EXPORT
CMR_ERROR CMRseymourCapture(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Releases a decomposition node, freeing it if this was the last reference.
 *
 * Decreases the reference counter by 1. If it reaches zero then it is freed. In that case, it is also called
 * recursively for the child nodes. Sets \p *pnode to \c NULL to prevent accidental usage.
 */

CMR_EXPORT
CMR_ERROR CMRseymourRelease(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE** pnode  /**< Pointer to Seymour decomposition node. \p *pnode is set to \c NULL. */
);

/**
 * \brief Creates an unknown decomposition node as a root.
 *
 * Copies \p matrix into the node.
 */

CMR_EXPORT
CMR_ERROR CMRseymourCreate(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE** pnode, /**< Pointer for storing the Seymour decomposition node. */
  bool isTernary,           /**< Whether we consider ternary matrices. */
  CMR_CHRMAT* matrix        /**< The matrix corresponding to this node; will be copied. */
);

/**
 * \brief Clones a decomposition node \p node into \p *pclone which represents the same matrix but has type
 *        \ref CMR_SEYMOUR_NODE_TYPE_UNKNOWN type and no child nodes.
 */

CMR_EXPORT
CMR_ERROR CMRseymourCloneUnknown(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node,   /**< Seymour decomposition node. */
  CMR_SEYMOUR_NODE** pclone /**< Pointer for storing the clone. */
);

/**
 * \brief Clones the union of subtrees of a Seymour decomposition, returning the copies.
 *
 * The set of Seymour decomposition nodes that are (grand-)children of any of the \p numSubtrees nodes \p subtreeRoots
 * is cloned. The respective clones are returned in \p clonedSubtrees.
 */

CMR_EXPORT
CMR_ERROR CMRseymourCloneSubtrees(CMR* cmr, size_t numSubtrees, CMR_SEYMOUR_NODE** subtreeRoots,
  CMR_SEYMOUR_NODE** clonedSubtrees);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_SEYMOUR_H */

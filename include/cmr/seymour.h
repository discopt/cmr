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
#include <cmr/network.h>
#include <cmr/series_parallel.h>

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup Seymour Seymour decomposition
 *
 * @{
 */


/**
 * \brief Flags that indicate how to decompose as a \f$ 3 \f$-sum.
 *
 * \see The desired types can be set by modifying \ref CMR_SEYMOUR_PARAMS.decomposeStrategy.
 **/

typedef enum
{
  CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_MASK = 15,
    /**< Bitmask for distributed rank treatment. */
  CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_PIVOT = 1,
    /**< Indicate to pivot a distributed rank distribution to a concentrated one; \see \ref seymour_decomposition. */
  CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_DELTASUM = 2,
    /**< Indicate to carry out a \f$ \Delta \f$-sum for distributed ranks; \see \ref seymour_decomposition. */
  CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM = 3,
    /**< Indicate to carry out a Y-sum for distributed ranks; \see \ref seymour_decomposition. */
  CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_MASK = 240,
    /**< Bitmask for distributed rank treatment. */
  CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT = 16,
    /**< Indicate to pivot a concentrated rank distribution to a distributed one; \see \ref seymour_decomposition. */
  CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM = 32,
    /**< Indicate to carry out a 3-sum for concentrated ranks; \see \ref seymour_decomposition. */
  CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_DELTASUM
    | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT,
    /**< This triggers only \f$ \Delta \f$-sum decompositions as defined by Seymour; \see \ref seymour_decomposition. */
  CMR_SEYMOUR_DECOMPOSE_FLAG_TRUEMPER = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_PIVOT
    | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM,
    /**< This triggers only 3-sum decompositions as defined by Truemper; \see \ref seymour_decomposition. */
} CMR_SEYMOUR_DECOMPOSE_FLAG;

/**
 * \brief Parameters for Seymour decomposition algorithm.
 */

typedef struct
{
  bool stopWhenIrregular;
  /**< \brief Whether to stop decomposing once irregularity is determined. */
  bool stopWhenNongraphic;
  /**< \brief Whether to stop decomposing once non-graphicness (or being non-network) is determined. */
  bool stopWhenNoncographic;
  /**< \brief Whether to stop decomposing once non-cographicness (or being non-conetwork) is determined. */
  bool stopWhenNeitherGraphicNorCoGraphic;
  /**< \brief Whether to stop decomposing once non-graphicness and non-cographicness (or not being network and not
   *          being conetwork) is determined. */

  bool seriesParallel;
  /**< \brief Whether to allow series-parallel operations in the decomposition tree; default: \c true */
  bool planarityCheck;
  /**< \brief Whether minors identified as graphic should still be checked for cographicness; default: \c false. */
  bool directGraphicness;
  /**< \brief Whether to use fast graphicness routines; default: \c true */
  bool preferGraphicness;
  /**< \brief Whether to first test for (co)graphicness (or being (co)network) before applying series-parallel
   *          reductions. */
  int decomposeStrategy;
  /**< \brief How to deal with 3-separations.
   **
   ** The value is a bit-wise OR of two decisions, one per **rank distribution**:
   **   - \ref CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_MASK indicates what to do if ranks are 1 and 1.
   **   - \ref CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_MASK indicates what to do if ranks are 2 and 0.
   **
   ** The possible choices for **distributed ranks** (1 and 1) are:
   **   - \ref CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_PIVOT pivot such that the rank distribution becomes concentrated.
   **   - \ref CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_DELTASUM for the \f$ \Delta \f$-sum (default).
   **   - \ref CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM for the Y-sum.
   **
   ** The possible choices for **concentrated ranks** (2 and 0) are:
   **   - \ref CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT pivot such that the rank distribution becomes distributed.
   **   - \ref CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM for the 3-sum (default).
   **
   ** \see \ref seymour_decomposition for a description of these sums.
   **
   ** A decomposition as described by Seymour can be selected via \ref CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR.
   ** A decomposition as used by Truemper can be selected via \ref CMR_SEYMOUR_DECOMPOSE_FLAG_TRUEMPER.
   ** The default is to not carry out any pivots and choose Seymour's or Truemper's definition depending on the rank
   ** distribution. */

  bool constructLeafGraphs;
  /**< \brief Whether to construct (co)graphs for all leaf nodes that are (co)graphic or (co)network. */
  bool constructAllGraphs;
  /**< \brief Whether to construct (co)graphs for all nodes that are (co)graphic or (co)network. */
} CMR_SEYMOUR_PARAMS;

/**
 * \brief Initializes the default parameters for regularity testing.
 *
 * These are selected for minimum running time.
 */

CMR_EXPORT
CMR_ERROR CMRseymourParamsInit(
  CMR_SEYMOUR_PARAMS* params  /**< Pointer to parameters. */
);

/**
 * \brief Statistics for Seymour decomposition algorithm.
 */

typedef struct
{
  uint32_t totalCount;                  /**< Total number of invocations. */
  double totalTime;                     /**< Total time of all invocations. */
  CMR_SP_STATISTICS seriesParallel;     /**< Statistics for series-parallel algorithm. */
  CMR_GRAPHIC_STATISTICS graphic;       /**< Statistics for direct (co)graphic checks. */
  CMR_NETWORK_STATISTICS network;       /**< Statistics for direct (co)network checks. */
  uint32_t sequenceExtensionCount;      /**< Number of extensions of sequences of nested minors. */
  double sequenceExtensionTime;         /**< Time of extensions of sequences of nested minors. */
  uint32_t sequenceGraphicCount;        /**< Number (co)graphicness tests applied to sequence of nested minors. */
  double sequenceGraphicTime;           /**< Time of (co)graphicness tests applied to sequence of nested minors. */
  uint32_t enumerationCount;            /**< Number of calls to enumeration algorithm for candidate 3-separations. */
  double enumerationTime;               /**< Time of enumeration of candidate 3-separations. */
  uint32_t enumerationCandidatesCount;  /**< Number of enumerated candidates for 3-separations. */
} CMR_SEYMOUR_STATS;


/**
 * \brief Initializes all statistics for Seymour decomposition computations.
 */

CMR_EXPORT
CMR_ERROR CMRseymourStatsInit(
  CMR_SEYMOUR_STATS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for Seymour decomposition computations.
 */

CMR_EXPORT
CMR_ERROR CMRseymourStatsPrint(
  FILE* stream,             /**< File stream to print to. */
  CMR_SEYMOUR_STATS* stats, /**< Pointer to statistics. */
  const char* prefix        /**< Prefix string to prepend to each printed line (may be \c NULL). */
);



struct _CMR_SEYMOUR_NODE;

typedef struct _CMR_SEYMOUR_NODE CMR_SEYMOUR_NODE;

typedef enum
{
  CMR_SEYMOUR_NODE_TYPE_IRREGULAR = -1,
    /**< Node represents 3-connected irregular matrix. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_UNKNOWN = 0,
    /**< Type of node is not yet determined. */
  CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL = 1,
    /**< Node represents a series-parallel reduction; has one child node. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_PIVOTS = 2,
    /**< Node represents an application of pivots; has one child node. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_GRAPH = 3,
    /**< Node represents a graphic or network leaf. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_COGRAPH = 4,
    /**< Node represents a cographic or conetwork leaf. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_PLANAR = 5,
    /**< Node represents a graphic and cographic (network and conetwork) leaf.  \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_R10 = 6,
    /**< Node represents a representation matrix of \f$ R_{10} \f$.  \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_ONESUM = 7,
    /**< Node represents a \f$ 1 \f$-sum of matrices; has at least 2 child nodes. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_TWOSUM = 8,
    /**< Node represents a \f$ 2 \f$-sum of matrices; has two child nodes. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_DELTASUM = 9,
    /**< Node represents a \f$ \Delta \f$-sum of matrices; has two child nodes. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_THREESUM = 10,
    /**< Node represents a \f$ 3 \f$-sum of matrices; has two child nodes. \see \ref seymour_decomposition. */
  CMR_SEYMOUR_NODE_TYPE_YSUM = 11
    /**< Node represents a Y-sum of matrices; has two child nodes. \see \ref seymour_decomposition. */
} CMR_SEYMOUR_NODE_TYPE;

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
 * \brief Returns the array of special rows of child \p childIndex.
 */

CMR_EXPORT
size_t* CMRseymourChildSpecialRows(
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  size_t childIndex       /**< Index of child to consider. */
);

/**
 * \brief Returns the array of special columns of child \p childIndex.
 */

CMR_EXPORT
size_t* CMRseymourChildSpecialColumns(
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
 * \brief Returns node's reference counter.
 */

CMR_EXPORT
size_t CMRseymourGetUsed(
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
 */

CMR_EXPORT
CMR_ERROR CMRseymourCreate(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE** pnode, /**< Pointer for storing the Seymour decomposition node. */
  bool isTernary,           /**< Whether we consider ternary matrices. */
  size_t numRows,           /**< Number of rows of represented matrix. */
  size_t numColumns         /**< Number of columns of represented matrix. */
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

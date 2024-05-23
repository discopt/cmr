#ifndef CMR_MATROID_H
#define CMR_MATROID_H

/**
 * \file matroid.h
 *
 * \author Matthias Walter
 *
 * \brief Functionality for matroids.
 */

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup Matroid Matroid decomposition
 *
 * @{
 */

/**
 * \brief Applie a pivot to \p matrix and returns the resulting matrix in \p *presult.
 *
 * Calculations are done over the binary field.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatBinaryPivot(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Matrix to work with. */
  size_t pivotRow,      /**< Row of the pivot. */
  size_t pivotColumn,   /**< Column of the pivot. */
  CMR_CHRMAT** presult  /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief Applies a sequence of pivots to \p matrix and returns the resulting matrix in \p *presult.
 *
 * Calculations are done over the ternary field.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatTernaryPivot(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Matrix to work with. */
  size_t pivotRow,      /**< Row of the pivot. */
  size_t pivotColumn,   /**< Column of the pivot. */
  CMR_CHRMAT** presult  /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief Applies a sequence of pivots to \p matrix and returns the resulting matrix in \p *presult.
 *
 * Calculations are done over the binary field.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatBinaryPivots(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Matrix to work with. */
  size_t numPivots,       /**< Number of pivots to carry out. */
  size_t* pivotRows,      /**< Array with rows of the pivots. */
  size_t* pivotColumns,   /**< Array with columns of the pivots. */
  CMR_CHRMAT** presult    /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief Applies a sequence of pivots to \p matrix and returns the resulting matrix in \p *presult.
 *
 * Calculations are done over the ternary field.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatTernaryPivots(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Matrix to work with. */
  size_t numPivots,       /**< Number of pivots to carry out. */
  size_t* pivotRows,      /**< Array with rows of the pivots. */
  size_t* pivotColumns,   /**< Array with columns of the pivots. */
  CMR_CHRMAT** presult    /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief A minor of a matroid.
 *
 * Specified by a sequence of pivots and a submatrix.
 */

typedef struct
{
  size_t numPivots;               /**< Number of pivots to apply. */
  size_t* pivotRows;              /**< Array with pivot rows. */
  size_t* pivotColumns;           /**< Array with pivot columns. */
  CMR_SUBMAT* remainingSubmatrix; /**< Submatrix that one finally needs to look at. */
} CMR_MINOR;

/**
 * \brief Creates a minor, allocating space for \p numPivots pivots and a remaining \p submatrix.
 */

CMR_EXPORT
CMR_ERROR CMRminorCreate(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_MINOR** pminor,   /**< Pointer for storing the minor. */
  size_t numPivots,     /**< Number of pivots. */
  CMR_SUBMAT* submatrix /**< Submatrix (may be \c NULL; is not copied). */
);

/**
 * \brief Frees the minor \p *pminor (if \p pminor is not \c NULL).
 */

CMR_EXPORT
CMR_ERROR CMRminorFree(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_MINOR** pminor  /**< Pointer to the minor (may be \c NULL). */
);

/**
 * \brief Writes the minor \p minor to the file \fileName by means of lists of row and column indices as well as
 *        pivot entries.
 */

CMR_EXPORT
CMR_ERROR CMRminorPrint(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_MINOR* minor,   /**< Minor to write. */
  size_t numRows,     /**< Number of rows of original matrix. */
  size_t numColumns,  /**< Number of columns of original matrix. */
  FILE* stream        /**< File stream to save minor to.. */
);

/**
 * \brief Writes the minor \p minor to the file \fileName by means of lists of row and column indices as well as
 *        pivot entries.
 */

CMR_EXPORT
CMR_ERROR CMRminorWriteToFile(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_MINOR* minor,     /**< Minor to write. */
  size_t numRows,       /**< Number of rows of original matrix. */
  size_t numColumns,    /**< Number of columns of original matrix. */
  const char* fileName  /**< File name to save minor to; \c NULL indicates stdout. */
);


struct _CMR_MATROID_DEC;

typedef struct _CMR_MATROID_DEC CMR_MATROID_DEC;

typedef enum
{
  CMR_MATROID_DEC_TYPE_IRREGULAR = -1,
    /**< Node represents 3-connected irregular minor. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_UNKNOWN = 0,
    /**< Type of node is not yet determined. */
  CMR_MATROID_DEC_TYPE_ONE_SUM = 1,
    /**< Node represents a \f$ 1 \f$-sum of matrices with an arbitrary number of child nodes.
     ** \see \ref matroid_decomposition.
     **/
  CMR_MATROID_DEC_TYPE_TWO_SUM = 2,
    /**< Node represents a 2-sum of matrices; has two child nodes. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_THREE_SUM = 3,
    /**< Node represents a 3-sum of matrices; has two child nodes. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_SERIES_PARALLEL = 4,
    /**< Node represents a series-parallel reduction; one child node. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_PIVOTS = 5,
    /**< Node represents an application of pivots; one child node. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_SUBMATRIX = 6,
    /**< Node represents the consideration of a submatrix; one child node. \see \ref matroid_decomposition. */

  CMR_MATROID_DEC_TYPE_GRAPH = 7,
    /**< Node represents a graphic leaf minor; one optional child node for non-cographic minor. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_COGRAPH = 8,
    /**< Node represents a cographic leaf minor; one optional child node for non-graphic minor. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_PLANAR = 9,
    /**< Node represents a planar (graphic and cographic) leaf minor; no child nodes. \see \ref matroid_decomposition. */

  CMR_MATROID_DEC_TYPE_R10 = -2,
    /**< Node represents a representation matrix of \f$ R_{10} \f$. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_FANO = -3,
    /**< Node represents a representation matrix of \f$ F_7 \f$. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_FANO_DUAL = -4,
    /**< Node represents a representation matrix of \f$ F_7^\star \f$. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_K5 = -5,
    /**< Node represents a representation matrix of \f$ M(K_5) \f$. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_K5_DUAL = -6,
    /**< Node represents a representation matrix of \f$ M(K_5)^\star \f$. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_K33 = -7,
    /**< Node represents a representation matrix of \f$ M(K_{3,3}) \f$. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_K33_DUAL = -8,
    /**< Node represents a representation matrix of \f$ M(K_{3,3})^\star \f$. \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_TYPE_DETERMINANT = -9,
    /**< Node represents a square matrix \f$ M \f$ with \f$ |\det(M)| = 2 \f$. \see \ref matroid_decomposition. */
} CMR_MATROID_DEC_TYPE;

/**
 * \brief Flags that indicate the type of \f$ 3 \f$-separation.
 *
 * \see The desired types can be set by modifying \ref CMR_REGULAR_PARAMS.threeSumStrategy.
 **/

typedef enum
{
  CMR_MATROID_DEC_THREESUM_FLAG_NO_PIVOTS = 0,
    /**< Indicate to not change the rank distribution; only serves as an option; each constructed node will have either
     **  \ref CMR_MATROID_DEC_THREESUM_FLAG_DISTRIBUTED_RANKS or
     **  \ref CMR_MATROID_DEC_THREESUM_FLAG_CONCENTRATED_RANK set. */
  CMR_MATROID_DEC_THREESUM_FLAG_DISTRIBUTED_RANKS = 1,
    /**< The two off-diagonal submatrices both have rank 1.
     **  \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_THREESUM_FLAG_CONCENTRATED_RANK = 2,
    /**< The bottom-left submatrix has rank 2 and the top-right submatrix has rank 0.
     **  \see \ref matroid_decomposition. */

  CMR_MATROID_DEC_THREESUM_FLAG_FIRST_WIDE = 4,
    /**< The first child node is of the form \f$ M_1^{\text{wide}} \f$; valid for distributed ranks.
     **  \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_THREESUM_FLAG_FIRST_TALL = 8,
    /**< The first child node is of the form \f$ M_1^{\text{tall}} \f$; valid for distributed ranks.
     **  \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_THREESUM_FLAG_FIRST_MIXED = 64,
    /**< The first child node is of the form \f$ M_1^{\text{mixed}} \f$; valid for concentrated rank.
     **  \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_THREESUM_FLAG_FIRST_ALLREPR = 128,
    /**< The first child node is of the form \f$ M_1^{\text{all-repr}} \f$; valid for concentrated rank.
     **  \see \ref matroid_decomposition. */

  CMR_MATROID_DEC_THREESUM_FLAG_SECOND_WIDE = 16,
    /**< The second child node is of the form \f$ M_2^{\text{wide}} \f$; valid for distributed ranks.
     **  \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_THREESUM_FLAG_SECOND_TALL = 32,
    /**< The second child node is of the form \f$ M_2^{\text{tall}} \f$; valid for distributed ranks.
     **  \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_THREESUM_FLAG_SECOND_MIXED = 256,
    /**< The second child node is of the form \f$ M_2^{\text{mixed}} \f$; valid for concentrated rank.
     **  \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_THREESUM_FLAG_SECOND_ALLREPR = 512,
    /**< The second child node is of the form \f$ M_2^{\text{all-repr}} \f$; valid for concentrated rank.
     **  \see \ref matroid_decomposition. */

  CMR_MATROID_DEC_THREESUM_FLAG_SEYMOUR = CMR_MATROID_DEC_THREESUM_FLAG_DISTRIBUTED_RANKS
    | CMR_MATROID_DEC_THREESUM_FLAG_FIRST_WIDE | CMR_MATROID_DEC_THREESUM_FLAG_SECOND_WIDE,
    /**< This combination of flags indicates a \f$ 3 \f$-sum as defined by Seymour.
     **  \see \ref matroid_decomposition. */
  CMR_MATROID_DEC_THREESUM_FLAG_TRUEMPER = CMR_MATROID_DEC_THREESUM_FLAG_CONCENTRATED_RANK
    | CMR_MATROID_DEC_THREESUM_FLAG_FIRST_MIXED | CMR_MATROID_DEC_THREESUM_FLAG_SECOND_MIXED,
    /**< This combination of flags indicates a \f$ 3 \f$-sum as defined by Truemper.
     **  \see \ref matroid_decomposition. */
} CMR_MATROID_DEC_THREESUM_FLAG;

/**
 * \brief Returns \c true iff the decomposition is over \f$ \mathbb{F}_3 \f$.
 */

CMR_EXPORT
bool CMRmatroiddecIsTernary(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns \c true iff the 3-sum decomposition node has \f$ 3 \f$-separation with two rank-1 matrices.
 */

CMR_EXPORT
bool CMRmatroiddecThreeSumDistributedRanks(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns \c true iff the 3-sum decomposition node has \f$ 3 \f$-separation with one rank-2 matrix.
 */

CMR_EXPORT
bool CMRmatroiddecThreeSumConcentratedRank(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns \c true iff the transposed matrix of the decomposition node \p dec is stored.
 */

CMR_EXPORT
bool CMRmatroiddecHasTranspose(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the matrix of the decomposition node \p dec (or \c NULL if it is not stored).
 */

CMR_EXPORT
CMR_CHRMAT* CMRmatroiddecGetMatrix(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the transposed matrix of the decomposition node \p dec (or \c NULL if it is not stored).
 */

CMR_EXPORT
CMR_CHRMAT* CMRmatroiddecGetTranspose(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the number of children of the decomposition node \p dec.
 */

CMR_EXPORT
size_t CMRmatroiddecNumChildren(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns a child of the decomposition node \p dec.
 */

CMR_EXPORT
CMR_MATROID_DEC* CMRmatroiddecChild(
  CMR_MATROID_DEC* dec,   /**< Decomposition node. */
  size_t childIndex       /**< Index of child. */
);

/**
 * \brief Returns the type of a decomposition node \p dec.
 */

CMR_EXPORT
CMR_MATROID_DEC_TYPE CMRmatroiddecType(
  CMR_MATROID_DEC* dec    /**< Decomposition node. */
);

/**
 * \brief Indicates graphicness/being network.
 *
 * Returns a positive value if the matrix corresponding to \p dec is graphic/network, zero if it is not known
 * and a negative value otherwise.
 */

CMR_EXPORT
int8_t CMRmatroiddecGraphicness(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Indicates cographicness/being conetwork.
 *
 * Returns a positive value if the matrix corresponding to \p dec is cographic/conetwork, zero if it is not known
 * and a negative value otherwise.
 */

CMR_EXPORT
int8_t CMRmatroiddecCographicness(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Indicates regularity/total unimodularity.
 *
 * Returns a positive value if the matrix corresponding to \p dec is regular/totally unimodular, zero if it is not known
 * and a negative value otherwise.
 */

CMR_EXPORT
int8_t CMRmatroiddecRegularity(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the number of rows.
 */

CMR_EXPORT
size_t CMRmatroiddecNumRows(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the number of columns.
 */

CMR_EXPORT
size_t CMRmatroiddecNumColumns(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the mapping of rows of child \p childIndex to this node's elements.
 */

CMR_EXPORT
CMR_ELEMENT* CMRmatroiddecChildRowsToParent(
  CMR_MATROID_DEC* dec, /**< Decomposition node. */
  size_t childIndex     /**< Index of child to consider. */
);

/**
 * \brief Returns the mapping of columns of child \p childIndex to this node's elements.
 */

CMR_EXPORT
CMR_ELEMENT* CMRmatroiddecChildColumnsToParent(
  CMR_MATROID_DEC* dec, /**< Decomposition node. */
  size_t childIndex     /**< Index of child to consider. */
);

/**
 * \brief Returns the graph (if available).
 */

CMR_EXPORT
CMR_GRAPH* CMRmatroiddecGraph(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the forest of the graph (if available).
 */

CMR_EXPORT
CMR_GRAPH_EDGE* CMRmatroiddecGraphForest(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the number of edges of the graph's forest (if available).
 */

CMR_EXPORT
size_t CMRmatroiddecGraphSizeForest(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the coforest of the graph (if available).
 */

CMR_EXPORT
CMR_GRAPH_EDGE* CMRmatroiddecGraphCoforest(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the number of edges of the graph's coforest (if available).
 */

CMR_EXPORT
size_t CMRmatroiddecGraphSizeCoforest(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns an array that indicates for the graph's edges whether they must be reversed (if available).
 */

CMR_EXPORT
bool* CMRmatroiddecGraphArcsReversed(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the cograph (if available).
 */

CMR_EXPORT
CMR_GRAPH* CMRmatroiddecCograph(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the number of edges of the cograph's forest (if available).
 */

CMR_EXPORT
size_t CMRmatroiddecCographSizeForest(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the forest of the cograph (if available).
 */

CMR_EXPORT
CMR_GRAPH_EDGE* CMRmatroiddecCographForest(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the number of edges of the cograph's coforest (if available).
 */

CMR_EXPORT
size_t CMRmatroiddecCographSizeCoforest(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the coforest of the cograph (if available).
 */

CMR_EXPORT
CMR_GRAPH_EDGE* CMRmatroiddecCographCoforest(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns an array that indicates for the cograph's edges whether they must be reversed (if available).
 */

CMR_EXPORT
bool* CMRmatroiddecCographArcsReversed(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the number of pivots (if available).
 */

CMR_EXPORT
size_t CMRmatroiddecNumPivots(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the array with the pivot rows (if available).
 */

CMR_EXPORT
size_t* CMRmatroiddecPivotRows(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Returns the array with the pivot columns (if available).
 */

CMR_EXPORT
size_t* CMRmatroiddecPivotColumns(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

/**
 * \brief Prints the decomposition \p dec to \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRmatroiddecPrint(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,     /**< Decomposition node. */
  FILE* stream,             /**< Stream to write to. */
  bool printChildren,       /**< Whether to recurse. */
  bool printParentElements, /**< Whether to print mapping of rows/columns to parent elements. */
  bool printMatrices,       /**< Whether to print matrices. */
  bool printGraphs,         /**< Whether to print graphs. */
  bool printReductions,     /**< Whether to print series-parallel reductions. */
  bool printPivots          /**< Whether to print pivots. */
);

/**
 * \brief Clones a decomposition node \p dec into \p *pclone which represents the same matrix but has type
 *        \ref CMR_MATROID_DEC_TYPE_UNKNOWN type.
 */

CMR_EXPORT
CMR_ERROR CMRmatroiddecCloneUnknown(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,     /**< The decomposition node. */
  CMR_MATROID_DEC** pclone  /**< Pointer for storing the clone. */
);

/**
 * \brief Increases the reference counter by 1.
 */

CMR_EXPORT
CMR_ERROR CMRmatroiddecCapture(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec  /**< Pointer to the decomposition node. */
);

/**
 * \brief Releases a decomposition node, freeing it if this was the last reference.
 *
 * Decreases the reference counter by 1. If it reaches zero then it is freed. In that case, it is also called
 * recursively for the child nodes.
 */

CMR_EXPORT
CMR_ERROR CMRmatroiddecRelease(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_MATROID_DEC** pdec  /**< Pointer to decomposition node. \p *pdec is set to \c NULL to prevent accidental usage. */
);

/**
 * \brief Creates an unknown decomposition node as a root.
 *
 * Copies \p matrix into the node.
 */

CMR_EXPORT
CMR_ERROR CMRmatroiddecCreateMatrixRoot(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_MATROID_DEC** pdec,   /**< Pointer for storing the decomposition node. */
  bool isTernary,           /**< Whether we consider ternary matrices. */
  CMR_CHRMAT* matrix        /**< The matrix corresponding to this node; will be copied. */
);

/**
 * \brief Clones the union of subtrees, returning the copies.
 *
 * The set of decomposition nodes that are (grand-)children of any of the \p numSubtrees nodes \p subtreeRoots is cloned.
 * The respective clones are returned in \p clonedSubtrees.
 */

CMR_EXPORT
CMR_ERROR CMRregularityCloneSubtrees(CMR* cmr, size_t numSubtrees, CMR_MATROID_DEC** subtreeRoots,
  CMR_MATROID_DEC** clonedSubtrees);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_MATROID_H */

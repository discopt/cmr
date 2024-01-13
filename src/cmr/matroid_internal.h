#ifndef CMR_MATROID_DEC_INTERNAL_H
#define CMR_MATROID_DEC_INTERNAL_H

#include <cmr/matroid.h>
#include <cmr/series_parallel.h>

struct _CMR_MATROID_DEC
{
  CMR_MATROID_DEC_TYPE type;                  /**< \brief Type of this node. */
  bool isTernary;                             /**< \brief Indicates whether this node belongs to a ternary matrix. */
  bool testedTwoConnected;                    /**< \brief Indicates that no 1-separation exists. */
  int8_t regularity;                          /**< \brief Matrix is (not) regular/totally unimodularity if positive
                                               **         (negative), or not determined if zero. */
  int8_t graphicness;                         /**< \brief Matrix is (not) graphic/network if positive
                                               **         (negative), or not determined if zero. */
  int8_t cographicness;                       /**< \brief Matrix is (not) cographic/conetwork if positive
                                               **         (negative), or not determined if zero. */
  bool testedR10;                             /**< \brief Matrix does not represent \f$ R_{10} \f$ unless \p type
                                               **         indicates this. */
  CMR_MATROID_DEC_THREESUM_TYPE threesumType; /**< \brief Type of 3-sum; an AND with
                                               **         \c CMR_MATROID_DEC_THREESUM_TYPE_DISTRIBUTED_RANKS
                                               **         indicates that the 3-sum has distributed ranks. */
  CMR_CHRMAT* matrix;                         /**< \brief Matrix representing this node. */
  CMR_CHRMAT* transpose;                      /**< \brief Tranpose of \ref matrix representing this node. */
  struct _CMR_MATROID_DEC* parent;            /**< \brief Parent node (\c NULL for decomposition root). */
  size_t numChildren;                         /**< \brief Number of child nodes. */
  struct _CMR_MATROID_DEC** children;         /**< \brief Array of child nodes. */

  size_t numRows;                             /**< \brief Length of \ref rowsParent. */
  size_t* rowsChild;                          /**< \brief Array for mapping each row to a row of the child (if
                                               **         applicable). */
  size_t* rowsParent;                         /**< \brief Array for mapping rows to rows of parent. */
  CMR_ELEMENT* rowsRootElement;               /**< \brief Array for mapping rows to elements of the root matrix. */

  size_t numColumns;                          /**< \brief Length of \ref columnsParent. */
  size_t* columnsChild;                       /**< \brief Array for mapping each column to a column of the child (if
                                               **        applicable). */
  size_t* columnsParent;                      /**< \brief Array for mapping columns to columns of parent. */
  CMR_ELEMENT* columnsRootElement;            /**< \brief Array for mapping rows to elements of the root matrix. */

  CMR_GRAPH* graph;                           /**< \brief Graph represented by this matrix. */
  CMR_GRAPH_EDGE* graphForest;                /**< \brief Array with edges of spanning forest of graph. */
  CMR_GRAPH_EDGE* graphCoforest;              /**< \brief Array with edges of coforest of graph. */
  bool* graphArcsReversed;                    /**< \brief Array indicating which arcs of the graph are reversed. */

  CMR_GRAPH* cograph;                         /**< \brief Graph represented by the transpose of this matrix. */
  CMR_GRAPH_EDGE* cographForest;              /**< \brief Array with edges of spanning forest of cograph. */
  CMR_GRAPH_EDGE* cographCoforest;            /**< \brief Array with edges of coforest of cograph. */
  bool* cographArcsReversed;                  /**< \brief Array indicating which arcs of the cograph are reversed. */

  bool testedSeriesParallel;                  /**< \brief Already searched for series-parallel reductions. */
  CMR_SP_REDUCTION* seriesParallelReductions; /**< \brief Array of series-parallel reductions. */
  size_t numSeriesParallelReductions;         /**< \brief Length of \p reductions. */

  size_t* pivotRows;                          /**< \brief Rows of pivots. */
  size_t* pivotColumns;                       /**< \brief Columns of pivots. */
  size_t numPivots;                           /**< \brief Number of pivots. */

  CMR_CHRMAT* nestedMinorsMatrix;             /**< \brief Equivalent binary matrix that displays the sequence of nested
                                               **         minors. */
  size_t* nestedMinorsSequenceNumRows;        /**< \brief Number of rows of sequence of nested minors. */
  size_t* nestedMinorsSequenceNumColumns;     /**< \brief Number of columns of sequence of nested minors. */
  size_t nestedMinorsLength;                  /**< \brief Length of sequence of nested minors. */
  CMR_ELEMENT* nestedMinorsRowsOriginal;      /**< \brief Maps rows of \p nestedMinorsMatrix to elements of
                                               **        \p matrix. */
  CMR_ELEMENT* nestedMinorsColumnsOriginal;   /**< \brief Maps columns of \p nestedMinorsMatrix to elements of
                                               **        \p matrix. */
};

/**
 * \brief Checks a matroid decomposition node for consistency.
 *
 * Checks all requirements defined in \ref matroids.
 *
 * \returns \c NULL if consistent. Otherwise, an explanation string is returned, which must freed with \c free().
 *
 * \see \ref CMRconsistencyAssert() for checking the returned string and aborting in case of inconsistency.
 */

char* CMRmatroiddecConsistency(
  CMR_MATROID_DEC* dec, /**< Decomposition. */
  bool recurse          /**< Whether all (grand-)children shall be checked, too. */
);

/**
 * \brief Creates a decomposition node as a child.
 *
 * Copies \p matrix and \p transpose into the node.
 */

CMR_ERROR CMRmatroiddecCreateChildFromMatrices(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_MATROID_DEC* parent,  /**< Parent node. */
  size_t childIndex,        /**< Child index of parent. */
  CMR_CHRMAT* matrix,       /**< The matrix corresponding to this node. */
  CMR_CHRMAT* transpose,    /**< The transpose matrix corresponding to this node. */
  size_t* rowsParent,       /**< Array for mapping rows to rows of parent. */
  size_t* columnsParent     /**< Array for mapping columns to columns of parent. */
);

/**
 * \brief Initialize an existing unknown decomposition node as a submatrix node whose child is of given \p type.
 */

CMR_ERROR CMRmatroiddecUpdateSubmatrix(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,     /**< Decomposition node. */
  CMR_SUBMAT* submatrix,    /**< Submatrix. */
  CMR_MATROID_DEC_TYPE type /**< Type of submatrix node. */
);

/**
 * \brief Initialize an existing unknown decomposition node as a 2-separation node according to the given \p separation.
 *
 * The two child nodes will be of type \ref CMR_MATROID_DEC_TYPE_UNKNOWN.
 */

CMR_ERROR CMRmatroiddecUpdateTwoSum(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec, /**< Decomposition node. */
  CMR_SEPA* separation  /**< 2-separation. */
);

/**
 * \brief Set regularity and (co)graphicness attributes of a decomposition tree.
 */

CMR_ERROR CMRmatroiddecSetAttributes(
  CMR_MATROID_DEC* dec  /**< Decomposition node. */
);

#endif /* CMR_MATROID_DEC_INTERNAL_H */

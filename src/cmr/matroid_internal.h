#ifndef CMR_MATROID_DEC_INTERNAL_H
#define CMR_MATROID_DEC_INTERNAL_H

#include "densematrix.h"

#include <cmr/matroid.h>
#include <cmr/seymour.h>
#include <cmr/series_parallel.h>


#ifdef __cplusplus
extern "C" {
#endif


struct _CMR_SEYMOUR_NODE
{
  size_t used;                                /**< \brief Reference counter. */
  CMR_SEYMOUR_NODE_TYPE type;                 /**< \brief Type of this node. */
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
  CMR_MATROID_DEC_THREESUM_FLAG threesumFlags;/**< \brief Type of 3-sum. */
  CMR_CHRMAT* matrix;                         /**< \brief Matrix representing this node. */
  CMR_CHRMAT* transpose;                      /**< \brief Tranpose of \ref matrix representing this node. */
  size_t numChildren;                         /**< \brief Number of child nodes. */
  struct _CMR_SEYMOUR_NODE** children;        /**< \brief Array of child nodes. */
  CMR_ELEMENT** childRowsToParent;            /**< \brief Array for mapping a child index to array of child rows to
                                               **         elements of this node. */
  CMR_ELEMENT** childColumnsToParent;         /**< \brief Array for mapping a child index to array of child columns to
                                               **         elements of this node. */

  size_t numRows;                             /**< \brief Length of \ref rowsParent. */
  size_t* rowsToChild;                         /**< \brief Array for mapping each row to a row of the child (if
                                               **         applicable). */

  size_t numColumns;                          /**< \brief Length of \ref columnsParent. */
  size_t* columnsToChild;                     /**< \brief Array for mapping each column to a column of the child (if
                                               **         applicable). */

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

  DenseBinaryMatrix* denseMatrix;             /**< \brief Dense support matrix for nested minors search. Rows and
                                               **         columns are not permuted, but pivots are applied. */
  CMR_ELEMENT* denseRowsOriginal;             /**< \brief Maps rows of \p denseMatrix to elements of \p matrix. */
  CMR_ELEMENT* denseColumnsOriginal;          /**< \brief Maps columns of \p denseMatrix to elements of \p matrix. */
  size_t* nestedMinorsRowsDense;              /**< \brief Maps rows of nested minor sequence to rows of \p dense. */
  size_t* nestedMinorsColumnsDense;           /**< \brief Maps columns of nested minor sequence to columns of
                                               **         \p dense. */

  size_t nestedMinorsLength;                  /**< \brief Length of sequence of nested minors. */
  size_t* nestedMinorsSequenceNumRows;        /**< \brief Number of rows of sequence of nested minors. */
  size_t* nestedMinorsSequenceNumColumns;     /**< \brief Number of columns of sequence of nested minors. */

  CMR_CHRMAT* nestedMinorsMatrix;             /**< \brief Sparse support matrix that displays the sequence of nested
                                               **         minors. Rows and columns are permuted accordingly. */
  CMR_CHRMAT* nestedMinorsTranspose;          /**< \brief Transpose of \p nestedMinorsMatrix. */
  CMR_ELEMENT* nestedMinorsRowsOriginal;      /**< \brief Maps rows of \p nestedMinorsMatrix to elements of
                                               **         \p matrix. */
  CMR_ELEMENT* nestedMinorsColumnsOriginal;   /**< \brief Maps columns of \p nestedMinorsDense to elements of
                                               **         \p matrix. */

  size_t nestedMinorsLastGraphic;             /**< \brief Last minor in sequence of nested minors that is graphic. */
  size_t nestedMinorsLastCographic;           /**< \brief Last minor in sequence of nested minors that is cographic. */
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
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  bool recurse            /**< Whether all (grand-)children shall be checked, too. */
);

/**
 * \brief Prints the decomposition \p child to \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRmatroiddecPrintChild(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* child,  /**< Seymour decomposition child node. */
  CMR_SEYMOUR_NODE* parent, /**< Seymour decomposition parent node. */
  size_t childIndex,        /**< Index of \p child as a child of \p parent. */
  FILE* stream,             /**< Stream to write to. */
  size_t indent,            /**< Indentation of this node. */
  bool printChildren,       /**< Whether to recurse. */
  bool printParentElements, /**< Whether to print mapping of rows/columns to parent elements (if \p parent is not \c NULL). */
  bool printMatrices,       /**< Whether to print matrices. */
  bool printGraphs,         /**< Whether to print graphs. */
  bool printReductions,     /**< Whether to print series-parallel reductions. */
  bool printPivots          /**< Whether to print pivots. */
);

/**
 * \brief Sets the number of children and allocates memory accordingly.
 */

CMR_ERROR CMRmatroiddecSetNumChildren(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  size_t numChildren      /**< Number of children. */
);

/**
 * \brief Initialize an existing unknown decomposition node as a 1-sum with \p numChildren children.
 */

CMR_ERROR CMRmatroiddecUpdateOneSum(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  size_t numChildren      /**< Number of child nodes. */
);

/**
 * \brief Creates a decomposition node as a child.
 *
 * Copies \p matrix and \p transpose into the node.
 */

CMR_ERROR CMRmatroiddecCreateChildFromMatrices(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* parent,     /**< Seymour decomposition parent node. */
  size_t childIndex,            /**< Child index of parent. */
  CMR_CHRMAT* matrix,           /**< The matrix corresponding to this node. */
  CMR_CHRMAT* transpose,        /**< The transpose matrix corresponding to this node. */
  CMR_ELEMENT* rowsToParent,    /**< Array for mapping rows to elements of parent. */
  CMR_ELEMENT* columnsToParent  /**< Array for mapping columns to elements of parent. */
);

/**
 * \brief Initialize an existing unknown decomposition node as a submatrix node whose child is of given \p type.
 */

CMR_ERROR CMRmatroiddecUpdateSubmatrix(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node,     /**< Seymour decomposition node. */
  CMR_SUBMAT* submatrix,      /**< Submatrix. */
  CMR_SEYMOUR_NODE_TYPE type  /**< Type of Seymour decomposition node. */
);

/**
 * \brief Initialize an existing unknown decomposition node as a 2-separation node according to the given \p separation.
 *
 * The two child nodes will be of type \ref CMR_MATROID_DEC_TYPE_UNKNOWN.
 */

CMR_ERROR CMRmatroiddecUpdateTwoSum(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  CMR_SEPA* separation    /**< 2-separation. */
);

/**
 * \brief Initialize an existing unknown decomposition node as a pivot node according to the given arrays of pivots.
 *
 * The unique child node will be of type \ref CMR_MATROID_DEC_TYPE_UNKNOWN.
 */

CMR_ERROR CMRmatroiddecUpdatePivots(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node, /**< Seymour decomposition node. */
  size_t numPivots,       /**< Number of pivots. */
  size_t* pivotRows,      /**< Array with pivot rows. */
  size_t* pivotColumns,   /**< Array with pivot columns. */
  CMR_CHRMAT* matrix,     /**< New matrix. */
  CMR_CHRMAT* transpose   /**< Transpose of \p matrix. */
);

/**
 * \brief Initialize an existing unknown decomposition node as a 3-sum node.
 *
 * The two child nodes remain uninitialized.
 */

CMR_ERROR CMRmatroiddecUpdateThreeSumInit(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

/**
 * \brief Creates wide first child of an initialized 3-sum node.
 *
 * A nonzero \p extraEntry indicates the bottom-right nonzero entry.
 */

CMR_ERROR CMRmatroiddecUpdateThreeSumCreateWideFirstChild(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node,     /**< Seymour decomposition node (initialized with \ref CMRmatroiddecUpdateThreeSumInit). */
  CMR_SEPA* separation,       /**< Separation. */
  size_t* rowsToChild,        /**< Array mapping rows to child rows. */
  size_t* columnsToChild,     /**< Array mapping columns to child columns. */
  size_t numChildBaseRows,    /**< Number of base rows of this child. */
  size_t numChildBaseColumns, /**< Number of base rows of this child. */
  size_t extraRow,            /**< Index of the extra row. */
  size_t extraColumn1,        /**< Index of 1st extra column. */
  size_t extraColumn2,        /**< Index of 2nd extra column, parallel to \p extraColumn1; equality is allowed. */
  int8_t extraEntry           /**< Sign of the extra entry. */
);

/**
 * \brief Creates wide second child of an initialized 3-sum node.
 *
 * A nonzero \p extraEntry indicates the bottom-right nonzero entry.
 */

CMR_ERROR CMRmatroiddecUpdateThreeSumCreateWideSecondChild(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node,     /**< Seymour decomposition node (initialized with \ref CMRmatroiddecUpdateThreeSumInit). */
  CMR_SEPA* separation,       /**< Separation. */
  size_t* rowsToChild,        /**< Array mapping rows to child rows. */
  size_t* columnsToChild,     /**< Array mapping columns to child columns. */
  size_t numChildBaseRows,    /**< Number of base rows of this child. */
  size_t numChildBaseColumns, /**< Number of base rows of this child. */
  size_t extraRow,            /**< Index of the extra row. */
  size_t extraColumn1,        /**< Index of 1st extra column. */
  size_t extraColumn2,        /**< Index of 2nd extra column, parallel to \p extraColumn1; equality is allowed. */
  int8_t extraEntry           /**< Sign of the extra entry, if known. */
);

/**
 * \brief Creates mixed first child of an initialized 3-sum node.
 *
 * A nonzero \p extraEntry indicates the bottom-right nonzero entry.
 */

CMR_ERROR CMRmatroiddecUpdateThreeSumCreateMixedFirstChild(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node,     /**< Seymour decomposition node (initialized with \ref CMRmatroiddecUpdateThreeSumInit). */
  CMR_SEPA* separation,       /**< Separation. */
  size_t* rowsToChild,        /**< Array mapping rows to child rows. */
  size_t* columnsToChild,     /**< Array mapping columns to child columns. */
  size_t numChildBaseRows,    /**< Number of base rows of this child. */
  size_t numChildBaseColumns, /**< Number of base rows of this child. */
  size_t extraRow1,           /**< Index of the first extra row. */
  size_t extraRow2,           /**< Index of the second extra row. */
  int8_t extraEntry           /**< Sign of the extra entry. */
);

/**
 * \brief Creates mixed second child of an initialized 3-sum node.
 *
 * A nonzero \p extraEntry indicates the top-left nonzero entry.
 */

CMR_ERROR CMRmatroiddecUpdateThreeSumCreateMixedSecondChild(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node,     /**< Seymour decomposition node (initialized with \ref CMRmatroiddecUpdateThreeSumInit). */
  CMR_SEPA* separation,       /**< Separation. */
  size_t* rowsToChild,        /**< Array mapping rows to child rows. */
  size_t* columnsToChild,     /**< Array mapping columns to child columns. */
  size_t numChildBaseRows,    /**< Number of base rows of this child. */
  size_t numChildBaseColumns, /**< Number of base rows of this child. */
  size_t extraColumn1,        /**< Index of the first extra column. */
  size_t extraColumn2,        /**< Index of the second extra column. */
  int8_t extraEntry           /**< Sign of the extra entry. */
);

/**
 * \brief Set regularity and (co)graphicness attributes of a decomposition tree.
 */

CMR_ERROR CMRmatroiddecSetAttributes(
  CMR_SEYMOUR_NODE* node  /**< Seymour decomposition node. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_MATROID_DEC_INTERNAL_H */

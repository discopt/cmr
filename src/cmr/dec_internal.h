#ifndef CMR_DEC_INTERNAL_H
#define CMR_DEC_INTERNAL_H

#include <cmr/dec.h>
#include <cmr/matroid.h>
#include <cmr/series_parallel.h>

struct _CMR_DEC
{
  CMR_DEC_TYPE type;                        /**< \brief Type of this node. */
  CMR_DEC_FLAGS flags;                      /**< \brief Flags for this node. */
  CMR_CHRMAT* matrix;                       /**< \brief Matrix representing this node. */
  CMR_CHRMAT* transpose;                    /**< \brief Tranpose of \ref matrix representing this node. */
  struct _CMR_DEC* parent;                  /**< \brief Parent node (\c NULL for decomposition root). */
  size_t numChildren;                       /**< \brief Number of child nodes. */
  struct _CMR_DEC** children;               /**< \brief Array of child nodes. */

  size_t numRows;                           /**< \brief Length of \ref rowsParent. */
  size_t* rowsParent;                       /**< \brief Array for mapping rows to rows of parent. */

  size_t numColumns;                        /**< \brief Length of \ref columnsParent. */
  size_t* columnsParent;                    /**< \brief Array for mapping columns to columns of parent. */

  CMR_GRAPH* graph;                         /**< \brief Graph represented by this matrix. */
  CMR_GRAPH_EDGE* graphForest;              /**< \brief Spanning forest of graph. */
  CMR_GRAPH_EDGE* graphCoforest;            /**< \brief Coforest of graph. */
  bool* graphArcsReversed;                  /**< \brief Array indicating which arcs of the graph are reversed. */

  CMR_GRAPH* cograph;                       /**< \brief Graph represented by the transpose of this matrix. */
  CMR_GRAPH_EDGE* cographForest;            /**< \brief Spanning forest of cograph. */
  CMR_GRAPH_EDGE* cographCoforest;          /**< \brief Coforest of cograph. */
  bool* cographArcsReversed;                /**< \brief Array indicating which arcs of the cograph are reversed. */

  CMR_SP_REDUCTION* reductions;             /**< \brief Array of series-parallel reductions. */
  size_t numReductions;                     /**< \brief Length of \p reductions. */

  CMR_CHRMAT* nestedMinorsMatrix;           /**< \brief Equivalent binary matrix that displays the sequence of nested
                                             **         minors. */
  size_t* nestedMinorsSequenceNumRows;      /**< \brief Number of rows of sequence of nested minors. */
  size_t* nestedMinorsSequenceNumColumns;   /**< \brief Number of columns of sequence of nested minors. */
  size_t nestedMinorsLength;                /**< \brief Length of sequence of nested minors. */
  CMR_ELEMENT* nestedMinorsRowsOriginal;    /**< \brief Maps rows of \p nestedMinorsMatrix to elements of \p matrix. */
  CMR_ELEMENT* nestedMinorsColumnsOriginal; /**< \brief Maps columns of \p nestedMinorsMatrix to elements of 
                                             **         \p matrix. */
};

/**
 * \brief Creates a decomposition node.
 */

CMR_ERROR CMRdecCreate(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DEC* parent,        /**< Parent decomposition node. */
  size_t numRows,         /**< Length of \p rowElements. */
  size_t* rowsParent,     /**< Array for mapping rows to rows of parent. */
  size_t numColumns,      /**< Length of \p columnsParent. */
  size_t* columnsParent,  /**< Array for mapping columns to columns of parent. */
  CMR_DEC** pdec          /**< Pointer for storing the created decomposition node. */
);

/**
 * \brief Construct \ref matrix and \ref tranpose based on \ref rowsParent and \ref columnsParent.
 */

CMR_ERROR CMRdecInheritMatrices(
  CMR* cmr,     /**< \ref CMR environment. */
  CMR_DEC* node /**< Decomposition node. */
);

/**
 * \brief Sets the number of child nodes and allocates memory.
 */

CMR_ERROR CMRdecSetNumChildren(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_DEC* node,      /**< Decomposition node. */
  size_t numChildren  /**< Number of children. */
);

/**
 * \brief Traverses decomposition tree to decide if \p node is regular, graphic or cographic.
 */

CMR_ERROR CMRdecComputeRegularity(
  CMR_DEC* node /**< Decomposition node. */
);

/**
 * \brief Translate the rows/columns of \p minor to the parent of \p node.
 */

CMR_ERROR CMRdecTranslateMinorToParent(
  CMR_DEC* node,    /**< Decomposition node. */
  CMR_MINOR* minor  /**< Minor to translate. */
);

/**
 * \brief Turns the given node into one for the given separation.
 */

CMR_ERROR CMRdecApplySeparation(
  CMR* cmr,       /**< \ref CMR environment. */
  CMR_DEC* dec,   /**< Decomposition node. */
  CMR_SEPA* sepa  /**< Separation. */
);

/**
 * \brief Prints the sequence of nested 3-connected minors for the matrix of a decomposition node.
 */

CMR_ERROR CMRdecPrintSequenceNested3ConnectedMinors(
  CMR* cmr,     /**< \ref CMR environment. */
  CMR_DEC* dec, /**< Decomposition node. */
  FILE* stream  /**< Stream to print to. */
);

#endif /* CMR_DEC_INTERNAL_H */

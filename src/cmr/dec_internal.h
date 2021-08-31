#ifndef CMR_DEC_INTERNAL_H
#define CMR_DEC_INTERNAL_H

#include <cmr/dec.h>
#include <cmr/matroid.h>
#include <cmr/series_parallel.h>

struct _CMR_DEC
{
  CMR_DEC_TYPE type;            /**< \brief Type of this node. */
  CMR_DEC_FLAGS flags;          /**< \brief Flags for this node. */
  CMR_CHRMAT* matrix;           /**< \brief Matrix representing this node. */
  CMR_CHRMAT* transpose;        /**< \brief Tranpose of \ref matrix representing this node. */
  struct _CMR_DEC* parent;      /**< \brief Parent node (\c NULL for decomposition root). */
  size_t numChildren;           /**< \brief Number of child nodes. */
  struct _CMR_DEC** children;   /**< \brief Array of child nodes. */

  size_t numRows;               /**< \brief Length of \ref rowsParent. */
  size_t* rowsParent;           /**< \brief Array for mapping rows to rows of parent. */

  size_t numColumns;            /**< \brief Length of \ref columnsParent. */
  size_t* columnsParent;        /**< \brief Array for mapping columns to columns of parent. */

  CMR_GRAPH* graph;             /**< \brief Graph represented by this matrix. */
  CMR_ELEMENT* edgeElements;    /**< \brief Array for mapping edges of \ref graph to elements. */

  CMR_GRAPH* cograph;           /**< \brief Graph represented by the transpose of this matrix. */
  CMR_ELEMENT* coedgeElements;  /**< \brief Array for mapping edges of \ref cograph to elements. */

  CMR_SP_REDUCTION* reductions; /**< \brief Array of series-parallel reductions. */
  size_t numReductions;         /**< \brief Length of \p reductions. */
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

#endif /* CMR_DEC_INTERNAL_H */

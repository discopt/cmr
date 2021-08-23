#ifndef CMR_TU_DECOMPOSITION_INTERNAL_H
#define CMR_TU_DECOMPOSITION_INTERNAL_H

#include <cmr/decomposition.h>

struct _CMR_TU_DEC
{
  CMR_TU_DEC_TYPE type;           /**< \brief Type of this node. */
  CMR_TU_DEC_FLAGS flags;         /**< \brief Flags for this node. */
  CMR_CHRMAT* matrix;          /**< \brief Matrix representing this node. */
  CMR_CHRMAT* transpose;       /**< \brief Tranpose of \ref matrix representing this node. */
  size_t numRows;             /**< \brief Length of \ref rowElements. */
  CMR_ELEMENT* rowElements;       /**< \brief Array for mapping rows to elements. */
  size_t* rowsParent;         /**< \brief Array for mapping rows to rows of parent. */
  size_t numColumns;          /**< \brief Length of \ref columnElements. */
  CMR_ELEMENT* columnElements;    /**< \brief Array for mapping rows to elements. */
  size_t* columnsParent;      /**< \brief Array for mapping columns to columns of parent. */
  size_t numPivots;           /**< \brief Number of pivots. */
  CMR_GRAPH* graph;            /**< \brief Graph represented by this matrix. */
  CMR_ELEMENT* edgeElements;      /**< \brief Array for mapping edges of \ref graph to elements. */
  CMR_GRAPH* cograph;          /**< \brief Graph represented by the transpose of this matrix. */
  CMR_ELEMENT* coedgeElements;    /**< \brief Array for mapping edges of \ref cograph to elements. */
  struct _CMR_TU_DEC* parent;     /**< \brief Parent node (\c NULL for decomposition root). */
  size_t numChildren;         /**< \brief Number of child nodes. */
  struct _CMR_TU_DEC** children;  /**< \brief Array of child nodes. */
};

CMR_ERROR CMRtudecCreate(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_TU_DEC* parent,         /**< Parent node. */
  size_t numRows,         /**< Length of \p rowElements. */
  size_t* rowsParent,     /**< Array for mapping rows to rows of parent. */
  size_t numColumns,      /**< Length of \p columnsParent. */
  size_t* columnsParent,  /**< Array for mapping columns to columns of parent. */
  CMR_TU_DEC** pdec           /**< Pointer for storing the created decomposition node. */
);

CMR_ERROR CMRtudecInheritElements(
  CMR* cmr,       /**< \ref CMR environment. */
  CMR_TU_DEC* node  /**< Node. */
);

CMR_ERROR CMRtudecInheritMatrices(
  CMR* cmr,       /**< \ref CMR environment. */
  CMR_TU_DEC* node  /**< Node. */
);

CMR_ERROR CMRtudecSetNumChildren(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_TU_DEC* node,       /**< Node. */
  size_t numChildren  /**< Number of children. */
);

CMR_ERROR CMRtudecComputeRegularity(
  CMR_TU_DEC* node  /**< Node. */
);

#endif /* CMR_TU_DECOMPOSITION_INTERNAL_H */

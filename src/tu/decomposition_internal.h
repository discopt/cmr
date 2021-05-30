#ifndef TU_DECOMPOSITION_INTERNAL_H
#define TU_DECOMPOSITION_INTERNAL_H

#include <tu/decomposition.h>

struct _TU_DEC
{
  TU_DEC_TYPE type;           /**< \brief Type of this node. */
  TU_DEC_FLAGS flags;         /**< \brief Flags for this node. */
  TU_CHRMAT* matrix;          /**< \brief Matrix representing this node. */
  TU_CHRMAT* transpose;       /**< \brief Tranpose of \ref matrix representing this node. */
  size_t numRows;             /**< \brief Length of \ref rowElements. */
  TU_ELEMENT* rowElements;       /**< \brief Array for mapping rows to elements. */
  size_t* rowsParent;         /**< \brief Array for mapping rows to rows of parent. */
  size_t numColumns;          /**< \brief Length of \ref columnElements. */
  TU_ELEMENT* columnElements;    /**< \brief Array for mapping rows to elements. */
  size_t* columnsParent;      /**< \brief Array for mapping columns to columns of parent. */
  size_t numPivots;           /**< \brief Number of pivots. */
  TU_GRAPH* graph;            /**< \brief Graph represented by this matrix. */
  TU_ELEMENT* edgeElements;      /**< \brief Array for mapping edges of \ref graph to elements. */
  TU_GRAPH* cograph;          /**< \brief Graph represented by the transpose of this matrix. */
  TU_ELEMENT* coedgeElements;    /**< \brief Array for mapping edges of \ref cograph to elements. */
  struct _TU_DEC* parent;     /**< \brief Parent node (\c NULL for decomposition root). */
  size_t numChildren;         /**< \brief Number of child nodes. */
  struct _TU_DEC** children;  /**< \brief Array of child nodes. */
};

TU_ERROR TUdecCreate(
  TU* tu,                 /**< \ref TU environment. */
  TU_DEC* parent,         /**< Parent node. */
  size_t numRows,         /**< Length of \p rowElements. */
  size_t* rowsParent,     /**< Array for mapping rows to rows of parent. */
  size_t numColumns,      /**< Length of \p columnsParent. */
  size_t* columnsParent,  /**< Array for mapping columns to columns of parent. */
  TU_DEC** pdec           /**< Pointer for storing the created decomposition node. */
);

TU_ERROR TUdecInheritElements(
  TU* tu,       /**< \ref TU environment. */
  TU_DEC* node  /**< Node. */
);

TU_ERROR TUdecInheritMatrices(
  TU* tu,       /**< \ref TU environment. */
  TU_DEC* node  /**< Node. */
);

TU_ERROR TUdecSetNumChildren(
  TU* tu,             /**< \ref TU environment. */
  TU_DEC* node,       /**< Node. */
  size_t numChildren  /**< Number of children. */
);

TU_ERROR TUdecComputeRegularity(
  TU_DEC* node  /**< Node. */
);

#endif /* TU_DECOMPOSITION_INTERNAL_H */

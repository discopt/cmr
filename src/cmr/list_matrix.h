#ifndef CMR_LIST_MATRIX_INTERNAL_H
#define CMR_LIST_MATRIX_INTERNAL_H

#include <cmr/env.h>

#include "hashtable.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Nonzero of list matrix.
 */

typedef struct _ListMatrixNonzero
{
  struct _ListMatrixNonzero* left;  /**< \brief Pointer to previous nonzero in the same row. */
  struct _ListMatrixNonzero* right; /**< \brief Pointer to next nonzero in the same row. */
  struct _ListMatrixNonzero* above; /**< \brief Pointer to previous nonzero in the same column. */
  struct _ListMatrixNonzero* below; /**< \brief Pointer to next nonzero in the same column. */
  size_t row;                       /**< \brief Row. */
  size_t column;                    /**< \brief Column. */
  char value;                       /**< \brief Matrix entry. */
  long special : 56;                /**< \brief Remaining bits (on 64 bit) may be used for a special purpose. */
} ListMatrixNonzero;

typedef struct
{
  ListMatrixNonzero head;             /**< \brief Dummy nonzero in that row/column. */
  size_t numNonzeros;                 /**< \brief Number of nonzeros in that row/column. */
  long long hashValue;                /**< \brief Hash value of this element. */
  CMR_LISTHASHTABLE_ENTRY hashEntry;  /**< \brief Entry in row or column hashtable. */
  size_t distance;                    /**< \brief Distance in breadth-first search. */
  size_t predecessor;                 /**< \brief Index of predecessor element in breadth-first search. */
  bool inQueue : 1;                   /**< \brief Whether this element is in a queue (e.g., of breadth-first search). */
//   long special : 63; 
  
  
  char lastBFS;                       /**< \brief Last breadth-first search that found this node.
                                       **< Is 0 initially, positive for search runs, -1 if marked and -2 for SP-reduced
                                       **< element. */
  bool specialBFS;                    /**< \brief Whether this is a special node in breadth-first search. */
} ListMatrixElement;



#ifdef __cplusplus
}
#endif

#endif /* CMR_LIST_MATRIX_INTERNAL_H */

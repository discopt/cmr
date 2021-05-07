#ifndef TU_HEAP_INTERNAL_H
#define TU_HEAP_INTERNAL_H

#include <tu/env.h>
#include "env_internal.h"
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Structure for min-heap with unsigned int values.
 */

typedef struct
{
  int size;       /**< \brief Current size of the heap. */
  int memKeys;    /**< \brief Memory for keys. */
  int* values;    /**< \brief Array that maps keys to values. */
  int* positions; /**< \brief Array that maps keys to heap positions. */
  int* data;      /**< \brief Array that maps heap positions to keys. */
} TU_INTHEAP;

/**
 * \brief Initializes an empty heap using stack memory.
 */

TU_ERROR TUintheapInitStack(
  TU* tu,           /**< \ref TU environment. */
  TU_INTHEAP* heap, /**< Heap pointer. */
  int memKeys       /**< Maximum number of elements and bound on key entries. */
);

/**
 * \brief Clears the given \p heap.
 */

TU_ERROR TUintheapClearStack(
  TU* tu,           /**< \ref TU environment. */
  TU_INTHEAP* heap  /**< Heap pointer. */
);

/**
 * \brief Inserts a \p key \p value pair into the heap.
 */

TU_ERROR TUintheapInsert(
  TU_INTHEAP* heap, /**< Heap pointer. */
  int key,          /**< Key of new element. */
  int value         /**< Value of new element. */
);

/**
 * \brief Decreases the value of \p key to \p newValue.
 */

TU_ERROR TUintheapDecrease(
  TU_INTHEAP* heap, /**< Heap pointer. */
  int key,          /**< Key of element. */
  int newValue      /**< New value of element. */
);

/**
 * \brief Decreases the value of \p key to \p newValue or inserts it.
 */

TU_ERROR TUintheapDecreaseInsert(
  TU_INTHEAP* heap, /**< Heap pointer. */
  int key,          /**< Key of element. */
  int newValue      /**< New value of element. */
);

/**
 * \brief Returns \c true if the heap is empty.
 */

static inline
bool TUintheapEmpty(
  TU_INTHEAP* heap  /**< Heap pointer. */
)
{
  return heap->size == 0;
}

/**
 * \brief Returns the key of the minimum element.
 */

static inline
int TUintheapMinimumKey(
  TU_INTHEAP* heap  /**< Heap pointer. */
)
{
  return heap->data[0];
}

/**
 * \brief Returns the value of the minimum element.
 */

static inline
int TUintheapMinimumValue(
  TU_INTHEAP* heap  /**< Heap. */
)
{
  return heap->values[heap->data[0]];
}

/**
 * \brief Returns \c true if an element with \p key is present in the heap.
 */

static inline
bool TUintheapContains(
  TU_INTHEAP* heap, /**< Heap pointer. */
  int key           /**< Key whose existence shall be checked. */
)
{
  return heap->positions[key] >= 0;
}

static inline
int TUintheapGetValue(
  TU_INTHEAP* heap, /*< Heap. */
  int key           /*< Key whose value shall be returned. */
)
{
  return heap->values[key];
}

/**
 * \brief Reutrns the value of \p key or \c INT_MAX if there is no such element.
 */

static inline
int TUintheapGetValueInfinity(
  TU_INTHEAP* heap, /**< Heap pointer. */
  int key           /**< Key to be searched. */
)
{
  if (heap->positions[key] >= 0)
    return heap->values[key];
  else
    return INT_MAX;
}

/**
 * \brief Extracts the minimum element and returns its key.
 */

int TUintheapExtractMinimum(
  TU_INTHEAP* heap  /**< Heap pointer. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_HEAP_INTERNAL_H */

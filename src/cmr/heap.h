#ifndef CMR_HEAP_INTERNAL_H
#define CMR_HEAP_INTERNAL_H

#include <cmr/env.h>
#include "env_internal.h"
#include <limits.h>

#ifdef __cplus
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
} CMR_INTHEAP;

/**
 * \brief Initializes an empty heap using stack memory.
 */

CMR_ERROR CMRintheapInitStack(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_INTHEAP* heap, /**< Heap pointer. */
  int memKeys       /**< Maximum number of elements and bound on key entries. */
);

/**
 * \brief Clears the given \p heap.
 */

CMR_ERROR CMRintheapClearStack(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_INTHEAP* heap  /**< Heap pointer. */
);

/**
 * \brief Inserts a \p key \p value pair into the heap.
 */

CMR_ERROR CMRintheapInsert(
  CMR_INTHEAP* heap, /**< Heap pointer. */
  int key,          /**< Key of new element. */
  int value         /**< Value of new element. */
);

/**
 * \brief Decreases the value of \p key to \p newValue.
 */

CMR_ERROR CMRintheapDecrease(
  CMR_INTHEAP* heap, /**< Heap pointer. */
  int key,          /**< Key of element. */
  int newValue      /**< New value of element. */
);

/**
 * \brief Decreases the value of \p key to \p newValue or inserts it.
 */

CMR_ERROR CMRintheapDecreaseInsert(
  CMR_INTHEAP* heap, /**< Heap pointer. */
  int key,          /**< Key of element. */
  int newValue      /**< New value of element. */
);

/**
 * \brief Returns \c true if the heap is empty.
 */

static inline
bool CMRintheapEmpty(
  CMR_INTHEAP* heap  /**< Heap pointer. */
)
{
  return heap->size == 0;
}

/**
 * \brief Returns the key of the minimum element.
 */

static inline
int CMRintheapMinimumKey(
  CMR_INTHEAP* heap  /**< Heap pointer. */
)
{
  return heap->data[0];
}

/**
 * \brief Returns the value of the minimum element.
 */

static inline
int CMRintheapMinimumValue(
  CMR_INTHEAP* heap  /**< Heap. */
)
{
  return heap->values[heap->data[0]];
}

/**
 * \brief Returns \c true if an element with \p key is present in the heap.
 */

static inline
bool CMRintheapContains(
  CMR_INTHEAP* heap, /**< Heap pointer. */
  int key           /**< Key whose existence shall be checked. */
)
{
  return heap->positions[key] >= 0;
}

static inline
int CMRintheapGetValue(
  CMR_INTHEAP* heap, /*< Heap. */
  int key           /*< Key whose value shall be returned. */
)
{
  return heap->values[key];
}

/**
 * \brief Reutrns the value of \p key or \c INT_MAX if there is no such element.
 */

static inline
int CMRintheapGetValueInfinity(
  CMR_INTHEAP* heap, /**< Heap pointer. */
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

int CMRintheapExtractMinimum(
  CMR_INTHEAP* heap  /**< Heap pointer. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_HEAP_INTERNAL_H */

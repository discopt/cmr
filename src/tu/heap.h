#ifndef TU_HEAP_INTERNAL_H
#define TU_HEAP_INTERNAL_H

#include <tu/env.h>
#include "env_internal.h"
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Structure for min-heap with unsigned int values.
 */

typedef struct
{
  int size;       /*< Current size of the heap. */
  int memKeys;    /*< Memory for keys. */
  int* values;    /*< Array that maps keys to values. */
  int* positions; /*< Array that maps keys to heap positions. */
  int* data;      /*< Array that maps heap positions to keys. */
} TU_INTHEAP;

void TUintheapInitStack(
  TU* tu,           /*< TU environment. */
  TU_INTHEAP* heap, /*< Heap. */
  int memKeys       /*< Memory for keys and heap. */
);

void TUintheapClearStack(
  TU* tu,           /*< TU environment. */
  TU_INTHEAP* heap  /*< Heap. */
);

void TUintheapInsert(
  TU_INTHEAP* heap, /*< Heap. */
  int key,          /*< Key of new element. */
  int value         /*< Value of new element. */
);

void TUintheapDecrease(
  TU_INTHEAP* heap, /*< Heap. */
  int key,          /*< Key of element. */
  int newValue      /*< New value of element. */
);

void TUintheapDecreaseInsert(
  TU_INTHEAP* heap, /*< Heap. */
  int key,          /*< Key of element. */
  int newValue      /*< New value of element. */
);

static inline
bool TUintheapEmpty(
  TU_INTHEAP* heap  /*< Heap. */
)
{
  return heap->size == 0;
}

static inline
int TUintheapMinimumKey(
  TU_INTHEAP* heap  /*< Heap. */
)
{
  return heap->data[0];
}

static inline
int TUintheapMinimumValue(
  TU_INTHEAP* heap  /*< Heap. */
)
{
  return heap->values[heap->data[0]];
}

static inline
bool TUintheapContains(
  TU_INTHEAP* heap, /*< Heap. */
  int key           /*< Some key. */
)
{
  return heap->positions[key] >= 0;
}

static inline
int TUintheapGetValue(
  TU_INTHEAP* heap, /*< Heap. */
  int key           /*< Some key. */
)
{
  return heap->values[key];
}

static inline
int TUintheapGetValueInfinity(
  TU_INTHEAP* heap, /*< Heap. */
  int key           /*< Some key. */
)
{
  if (heap->positions[key] >= 0)
    return heap->values[key];
  else
    return INT_MAX;
}

int TUintheapExtractMinimum(
  TU_INTHEAP* heap  /*< Heap. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_HEAP_INTERNAL_H */

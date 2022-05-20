// #define CMR_DEBUG /* Uncomment to debug the heap. */
// #define CMR_DEBUG_HEAP_CONTENT /* Uncomment to print the whole heap after each operation. */

#include "heap.h"

#include <assert.h>
#include <limits.h>

CMR_ERROR CMRintheapInitStack(CMR* cmr, CMR_INTHEAP* heap, int memKeys)
{
  assert(cmr);
  assert(heap);
  assert(memKeys > 0);

  heap->memKeys = memKeys;
  heap->size = 0;
  heap->positions = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &heap->positions, memKeys) );
  for (int i = 0; i < memKeys; ++i)
    heap->positions[i] = -1;
  heap->values = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &heap->values, memKeys) );
  heap->data = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &heap->data, memKeys) );

  return CMR_OKAY;
}

CMR_ERROR CMRintheapClearStack(CMR* cmr, CMR_INTHEAP* heap)
{
  assert(cmr);
  assert(heap);

  CMR_CALL( CMRfreeStackArray(cmr, &heap->data) );
  CMR_CALL( CMRfreeStackArray(cmr, &heap->values) );
  CMR_CALL( CMRfreeStackArray(cmr, &heap->positions) );
  heap->memKeys = 0;

  return CMR_OKAY;
}

#if defined(CMR_DEBUG_HEAP_CONTENT)
static
void debugHeap(CMR_INTHEAP* heap)
{
  printf("                    Heap:");
  for (int i = 0; i < heap->size; ++i)
  {
    printf(" %d:%d->%d", i, heap->data[i], heap->values[heap->data[i]]);
  }
  printf("\n");
  fflush(stdout);
}
#else
static inline
void debugHeap(CMR_INTHEAP* heap)
{
  CMR_UNUSED(heap);
}
#endif /* CMR_DEBUG_HEAP_CONTENT */

CMR_ERROR CMRintheapInsert(CMR_INTHEAP* heap, int key, int value)
{
  assert(heap);
  assert(key >= 0);
  assert(key < heap->memKeys);
  assert(heap->size < heap->memKeys);

  CMRdbgMsg(20, "Heap insert: %d->%d.\n", key, value);

  heap->data[heap->size] = key;
  heap->positions[key] = heap->size;
  heap->values[key] = value;

  int currentKey = key;
  int current = heap->size;
  int currentValue = value;
  while (current > 0)
  {
    int parent = (current-1) / 2;
    int parentKey = heap->data[parent];
    int parentValue = heap->values[parentKey];
    if (parentValue <= currentValue)
      break;

    /* Move current upwards. */
    heap->positions[currentKey] = parent;
    heap->positions[parentKey] = current;
    heap->data[parent] = currentKey;
    heap->data[current]  = parentKey;
    current = parent;
    currentKey = parentKey;
    currentValue = parentValue;
  }

  ++heap->size;

  debugHeap(heap);

  return CMR_OKAY;
}

CMR_ERROR CMRintheapDecrease(CMR_INTHEAP* heap, int key, int newValue)
{
  assert(heap);
  assert(heap->positions[key] >= 0);

  int currentKey = key;
  int current = heap->positions[currentKey];
  CMRdbgMsg(20, "Heap decrease: %d->%d to %d->%d.\n", key, heap->values[currentKey], key, newValue);
  heap->values[currentKey] = newValue;
  int currentValue = newValue;
  while (current > 0)
  {
    int parent = (current-1) / 2;
    int parentKey = heap->data[parent];
    int parentValue = heap->values[parentKey];
    CMRdbgMsg(12, "Parent: %d->%d.\n", parentKey, parentValue);
    if (parentValue <= currentValue)
      break;

    /* Move current upwards. */
    heap->positions[currentKey] = parent;
    heap->positions[parentKey] = current;
    heap->data[parent] = currentKey;
    heap->data[current]  = parentKey;
    current = parent;
    currentKey = parentKey;
    currentValue = parentValue;
  }

  debugHeap(heap);

  return CMR_OKAY;
}

CMR_ERROR CMRintheapDecreaseInsert(CMR_INTHEAP* heap, int key, int newValue)
{
  assert(heap);
  assert(key >= 0);
  assert(key < heap->memKeys);

  CMRdbgMsg(20, "Heap decrease-insert: %d->%d.\n", key, newValue);

  int currentKey = key;
  int current;
  if (heap->positions[key] >= 0)
  {
    current = heap->positions[currentKey];
  }
  else
  {
    current = heap->size;
    heap->size++;
    heap->positions[key] = current;
    heap->data[current] = key;
  }
  heap->values[currentKey] = newValue;
  int currentValue = newValue;
  CMRdbgMsg(22, "Initial is %d:%d->%d.\n", current, currentKey, currentValue);

  while (current > 0)
  {
    int parent = (current-1) / 2;
    int parentKey = heap->data[parent];
    int parentValue = heap->values[parentKey];
    CMRdbgMsg(22, "Parent: %d:%d->%d, child value: %d.\n", parent, parentKey, parentValue, currentValue);
    if (parentValue <= currentValue)
      break;

    /* Move current upwards. */
    heap->positions[currentKey] = parent;
    heap->positions[parentKey] = current;
    heap->data[parent] = currentKey;
    heap->data[current]  = parentKey;
    current = parent;
  }

  debugHeap(heap);

  return CMR_OKAY;
}

int CMRintheapExtractMinimum(CMR_INTHEAP* heap)
{
  assert(heap);

  int extracted = heap->data[0];
  heap->positions[extracted] = -1;
  heap->data[0] = heap->data[heap->size - 1];
  heap->positions[heap->data[0]] = 0;
  --heap->size;

  CMRdbgMsg(20, "Heap extract: %d->%d.\n", extracted, heap->values[extracted]);

  int current = 0;
  int currentKey = heap->data[0];
  int currentValue = heap->values[currentKey];
  while (2*current+1 < heap->size)
  {
    int left = 2*current + 1;
    int right = left+1;
    int leftKey = heap->data[left];
    int leftValue = heap->values[leftKey];
    int rightKey = right < heap->size ? heap->data[right] : -1;
    int rightValue = right < heap->size ? heap->values[rightKey] : INT_MAX;
    if (leftValue < rightValue)
    {
      if (currentValue <= leftValue)
        break;

      CMRdbgMsg(22, "Swapping %d:%d->%d with left child %d:%d->%d.\n", current, currentKey, currentValue, left, leftKey,
        leftValue);
      heap->positions[leftKey] = current;
      heap->positions[currentKey] = left;
      heap->data[current] = leftKey;
      heap->data[left] = currentKey;
      current = left;
    }
    else
    {
      if (currentValue <= rightValue)
        break;

      CMRdbgMsg(22, "Swapping %d:%d->%d with right child %d:%d->%d.\n", current, currentKey, currentValue, right, rightKey,
        rightValue);

      heap->positions[rightKey] = current;
      heap->positions[currentKey] = right;
      heap->data[current] = rightKey;
      heap->data[right] = currentKey;
      current = right;
    }
  }

  debugHeap(heap);

  return extracted;
}

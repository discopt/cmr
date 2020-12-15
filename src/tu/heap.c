#include "heap.h"

#include <assert.h>
#include <limits.h>

void TUintheapInitStack(TU* tu, TU_INTHEAP* heap, int memKeys)
{
  assert(tu);
  assert(heap);
  assert(memKeys > 0);

  heap->memKeys = memKeys;
  heap->size = 0;
  heap->positions = NULL;
  TUallocStackArray(tu, &heap->positions, memKeys);
  for (int i = 0; i < memKeys; ++i)
    heap->positions[i] = -1;
  heap->values = NULL;
  TUallocStackArray(tu, &heap->values, memKeys);
  heap->data = NULL;
  TUallocStackArray(tu, &heap->data, memKeys);
}

void TUintheapClearStack(TU* tu, TU_INTHEAP* heap)
{
  assert(tu);
  assert(heap);

  TUfreeStackArray(tu, &heap->data);
  TUfreeStackArray(tu, &heap->values);
  TUfreeStackArray(tu, &heap->positions);
  heap->memKeys = 0;
}

void TUintheapInsert(TU_INTHEAP* heap, int key, int value)
{
  assert(heap);
  assert(key >= 0);
  assert(key < heap->memKeys);

  assert(heap->size < heap->memKeys);
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
}

void TUintheapDecrease(TU_INTHEAP* heap, int key, int newValue)
{
  assert(heap);
  assert(heap->positions[key] >= 0);

  int currentKey = key;
  int current = heap->positions[currentKey];
  heap->values[currentKey] = newValue;
  int currentValue = newValue;
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
}

void TUintheapDecreaseInsert(TU_INTHEAP* heap, int key, int newValue)
{
  assert(heap);
  assert(key >= 0);
  assert(key < heap->memKeys);

  int currentKey = key;
  int current;
  if (heap->positions[key] >= 0)
    current = heap->positions[currentKey];
  else
  {
    current = heap->size;
    heap->size++;
    heap->positions[key] = current;
    heap->data[current] = key;
  }
  heap->values[currentKey] = newValue;
  int currentValue = newValue;
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
}

int TUintheapExtractMinimum(TU_INTHEAP* heap)
{
  assert(heap);

  int extracted = heap->data[0];
  heap->positions[extracted] = -1;
  heap->data[0] = heap->data[heap->size - 1];
  heap->positions[heap->data[0]] = 0;
  --heap->size;

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

      heap->positions[leftKey] = current;
      heap->positions[currentKey] = left;
      heap->data[current] = leftKey;
      heap->data[left] = currentKey;
      current = left;
      currentKey = leftKey;
      currentValue = leftValue;
    }
    else
    {
      if (currentValue <= rightValue)
        break;

      heap->positions[rightKey] = current;
      heap->positions[currentKey] = right;
      heap->data[current] = rightKey;
      heap->data[right] = currentKey;
      current = right;
      currentKey = rightKey;
      currentValue = rightValue;
    }
  }
  
  return extracted;
}

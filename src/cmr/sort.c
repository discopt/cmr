// #define TU_DEBUG /* Uncomment to debug the sorting. */

#include "sort.h"

#include <assert.h>
#include <stdlib.h>

#define MOVE_ELEMENT(source, target, size) \
  do \
  { \
    size_t __size = (size); \
    char* __s = (source); \
    char* __t = (target); \
    do \
    { \
      *__t++ = *__s++; \
    } \
    while (--__size > 0); \
  } \
  while (false)

CMR_ERROR TUsort(TU* tu, size_t length, void* array, size_t elementSize, int (*compare)(const void*, const void*))
{
  assert(tu);

  qsort(array, length, elementSize, compare);

  return CMR_OKAY;
}

CMR_ERROR TUsort2(TU* tu, size_t length, void* array1, size_t elementSize1, void* array2, size_t elementSize2,
  int (*compare)(const void**, const void**))
{
  assert(tu);

  TUassertStackConsistency(tu);

  /* Allocate space for temporary elements of both arrays. */
  char* temp1 = NULL;
  TU_CALL( TUallocStackArray(tu, &temp1, elementSize1) );
  char* temp2 = NULL;
  TU_CALL( TUallocStackArray(tu, &temp2, elementSize2) );

  /* Create an array with pointers into array1. */
  void** pointerArray = NULL;
  TU_CALL( TUallocStackArray(tu, &pointerArray, length) );
  void* pointer = array1;
  for (size_t i = 0; i < length; ++i)
  {
    pointerArray[i] = pointer;
    pointer += elementSize1;
  }

  qsort(pointerArray, length, sizeof(size_t*), (int (*)(const void*, const void*)) compare);

  /* Reorder array1 and array2 according to the array of pointers. */
  for (size_t i = 0; i < length; ++i)
  {
    TUdbgMsg(0, "i = %d, pointerArray[i] = %p array1 = %p, difference = %ld, scaled = %d\n", i, pointerArray[i], array1,
      pointerArray[i] - array1, (pointerArray[i] - array1) / elementSize1);
    if (i != (pointerArray[i] - array1) / elementSize1)
    {
      MOVE_ELEMENT(array1 + elementSize1 * i, temp1, elementSize1); /* temp1 = array1[i] */
      MOVE_ELEMENT(array2 + elementSize2 * i, temp2, elementSize2); /* temp2 = array2[i] */
      size_t j;
      size_t k = i;
      while (i != (j = (pointerArray[k] - array1) / elementSize1))
      {
        MOVE_ELEMENT(array1 + elementSize1 * j, array1 + elementSize1 * k, elementSize1); /* array1[k] = array1[j] */
        MOVE_ELEMENT(array2 + elementSize2 * j, array2 + elementSize2 * k, elementSize2); /* array2[k] = array2[j] */
        pointerArray[k] = array1 + elementSize1 * k;
        k = j;
      }
      MOVE_ELEMENT(temp1, array1 + elementSize1 * k, elementSize1); /* array1[k] = temp1 */
      MOVE_ELEMENT(temp2, array2 + elementSize2 * k, elementSize2); /* array2[k] = temp2 */
      pointerArray[k] = array1 + elementSize1 * k;
    }
  }

  /* Free the temporary space. */
  TU_CALL( TUfreeStackArray(tu, &pointerArray) );
  TU_CALL( TUfreeStackArray(tu, &temp2) );
  TU_CALL( TUfreeStackArray(tu, &temp1) );
  TUassertStackConsistency(tu);

  return CMR_OKAY;
}


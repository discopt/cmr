#define CMR_DEBUG /* Uncomment to debug the sorting. */

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

CMR_ERROR CMRsort(CMR* cmr, size_t length, void* array, size_t elementSize, int (*compare)(const void*, const void*))
{
  CMR_UNUSED(cmr);

  assert(cmr);

  qsort(array, length, elementSize, compare);

  return CMR_OKAY;
}

CMR_ERROR CMRsort2(CMR* cmr, size_t length, void* array1, size_t elementSize1, void* array2, size_t elementSize2,
  int (*compare)(const void**, const void**))
{
  assert(cmr);

  char* data1 = array1;
  char* data2 = array2;

  CMRassertStackConsistency(cmr);

  /* Allocate space for temporary elements of both arrays. */
  char* temp1 = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &temp1, elementSize1) );
  char* temp2 = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &temp2, elementSize2) );

  /* Create an array with pointers into array1. */
  char** pointerArray = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &pointerArray, length) );
  char* pointer = data1;
  for (size_t i = 0; i < length; ++i)
  {
    pointerArray[i] = pointer;
    pointer += elementSize1;
  }

  qsort(pointerArray, length, sizeof(size_t*), (int (*)(const void*, const void*)) compare);

  /* Reorder data1 and data2 according to the array of pointers. */
  for (size_t i = 0; i < length; ++i)
  {
    CMRdbgMsg(0, "i = %d, pointerArray[i] = %p array1 = %p, difference = %ld, scaled = %d\n", i, pointerArray[i], data1,
      pointerArray[i] - data1, (pointerArray[i] - data1) / elementSize1);
    if (i != (pointerArray[i] - data1) / elementSize1)
    {
      MOVE_ELEMENT(data1 + elementSize1 * i, temp1, elementSize1); /* temp1 = array1[i] */
      MOVE_ELEMENT(data2 + elementSize2 * i, temp2, elementSize2); /* temp2 = array2[i] */
      size_t j;
      size_t k = i;
      while (i != (j = (pointerArray[k] - data1) / elementSize1))
      {
        MOVE_ELEMENT(data1 + elementSize1 * j, data1 + elementSize1 * k, elementSize1); /* data1[k] = data1[j] */
        MOVE_ELEMENT(data2 + elementSize2 * j, data2 + elementSize2 * k, elementSize2); /* data2[k] = data2[j] */
        pointerArray[k] = data1 + elementSize1 * k;
        k = j;
      }
      MOVE_ELEMENT(temp1, data1 + elementSize1 * k, elementSize1); /* data1[k] = temp1 */
      MOVE_ELEMENT(temp2, data2 + elementSize2 * k, elementSize2); /* data2[k] = temp2 */
      pointerArray[k] = data1 + elementSize1 * k;
    }
  }

  /* Free the temporary space. */
  CMR_CALL( CMRfreeStackArray(cmr, &pointerArray) );
  CMR_CALL( CMRfreeStackArray(cmr, &temp2) );
  CMR_CALL( CMRfreeStackArray(cmr, &temp1) );
  CMRassertStackConsistency(cmr);

  return CMR_OKAY;
}


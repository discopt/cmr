#ifndef CMR_SORT_INTERNAL_H
#define CMR_SORT_INTERNAL_H

#include "env_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Sorts an array.
 *
 * The user must provide a \p compare function that is called several times with pointers into the first array.
 * The arguments must be dereferenced once in order to access the actual element.
 */

CMR_ERROR CMRsort(
  CMR* cmr,                                 /**< \ref CMR environment. */
  size_t length,                            /**< Number of elements in both arrays. */
  void *array,                              /**< Pointer to the first element of the array.  */
  size_t elementSize,                       /**< Size (in bytes) of each element of the array. */
  int (*compare)(const void*, const void*)  /**< Comparison function for comparing two elements. */
);

/**
 * \brief Sorts two arrays simultaneously.
 *
 * Sorts two arrays using the qsort function from the C library. The user must provide a \p compare function that is
 * called several times with pointers to pointers into the first array. Hence, the arguments must be dereferenced twice
 * in order to access the actual element.
 */

CMR_ERROR CMRsort2(
  CMR* cmr,                                     /**< \ref CMR environment. */
  size_t length,                              /**< Number of elements in both arrays. */
  void *array1,                               /**< Pointer to the first element of the first array.  */
  size_t elementSize1,                        /**< Size (in bytes) of each element of the first array. */
  void* array2,                               /**< Pointer to the first element of the second array.  */
  size_t elementSize2,                        /**< Size (in bytes) of each element of the second array. */
  int (*compare)(const void**, const void**)  /**< Comparison function for comparing two elements. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_SORT_INTERNAL_H */

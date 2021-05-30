#ifndef TU_ELEMENT_H
#define TU_ELEMENT_H

#include <tu/env.h>

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int TU_ELEMENT;

TU_EXPORT
const char* TUelementString(
  TU_ELEMENT element,  /**< Element to print. */
  char* buffer      /**< Buffer of size at least 32. May be \c NULL, in which case a static buffer is used. */
);

static inline
TU_ELEMENT TUrowToElement(
  size_t row  /**< Row index. */
)
{
  return -1 - (int)row;
}

static inline
TU_ELEMENT TUcolumnToElement(
  size_t column  /**< Column index. */
)
{
  return 1 + (int)column;
}

static inline
bool TUelementIsRow(
  TU_ELEMENT element /**< Element to check. */
)
{
  return element < 0;
}

static inline
size_t TUelementToRowIndex(
  TU_ELEMENT element /**< Element to convert. */
)
{
  assert(element < 0);
  return -1 - element;
}

static inline
bool TUelementIsColumn(
  TU_ELEMENT element /**< Element to check. */
)
{
  return element > 0;
}

static inline
size_t TUelementToColumnIndex(
  TU_ELEMENT element /**< Element to convert. */
)
{
  assert(element > 0);
  return -1 + element;
}

#ifdef __cplusplus
}
#endif

#endif /* TU_ELEMENT_H */

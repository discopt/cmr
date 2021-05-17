#ifndef TU_ELEMENT_H
#define TU_ELEMENT_H

#include <tu/env.h>

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int Element;

TU_EXPORT
const char* TUelementString(
  Element element,  /**< Element to print. */
  char* buffer      /**< Buffer of size at least 32. May be \c NULL, in which case a static buffer is used. */
);

static inline
Element TUrowToElement(
  size_t row  /**< Row index. */
)
{
  return -1 - (int)row;
}

static inline
Element TUcolumnToElement(
  size_t column  /**< Column index. */
)
{
  return 1 + (int)column;
}

static inline
bool TUelementIsRow(
  Element element /**< Element to check. */
)
{
  return element < 0;
}

static inline
size_t TUelementToRowIndex(
  Element element /**< Element to convert. */
)
{
  assert(element < 0);
  return -1 - element;
}

static inline
bool TUelementIsColumn(
  Element element /**< Element to check. */
)
{
  return element > 0;
}

static inline
size_t TUelementToColumnIndex(
  Element element /**< Element to convert. */
)
{
  assert(element > 0);
  return -1 + element;
}

#ifdef __cplusplus
}
#endif

#endif /* TU_ELEMENT_H */

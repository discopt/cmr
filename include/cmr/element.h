#ifndef CMR_ELEMENT_H
#define CMR_ELEMENT_H

#include <cmr/env.h>

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int CMR_ELEMENT;

CMR_EXPORT
const char* CMRelementString(
  CMR_ELEMENT element,  /**< Element to print. */
  char* buffer          /**< Buffer of size at least 32. May be \c NULL, in which case a static buffer is used. */
);

static inline
CMR_ELEMENT CMRrowToElement(
  size_t row  /**< Row index. */
)
{
  return -1 - (int)row;
}

static inline
CMR_ELEMENT CMRcolumnToElement(
  size_t column  /**< Column index. */
)
{
  return 1 + (int)column;
}

static inline
bool CMRelementIsRow(
  CMR_ELEMENT element /**< Element to check. */
)
{
  return element < 0;
}

static inline
size_t CMRelementToRowIndex(
  CMR_ELEMENT element /**< Element to convert. */
)
{
  assert(element < 0);
  return -1 - element;
}

static inline
bool CMRelementIsColumn(
  CMR_ELEMENT element /**< Element to check. */
)
{
  return element > 0;
}

static inline
size_t CMRelementToColumnIndex(
  CMR_ELEMENT element /**< Element to convert. */
)
{
  assert(element > 0);
  return -1 + element;
}

#ifdef __cplusplus
}
#endif

#endif /* CMR_ELEMENT_H */

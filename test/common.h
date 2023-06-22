#ifndef CMR_TEST_COMMON_H
#define CMR_TEST_COMMON_H

#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

CMR_ERROR stringToDoubleMatrix(CMR* cmr, CMR_DBLMAT** matrix, const char* string);

CMR_ERROR stringToIntMatrix(CMR* cmr, CMR_INTMAT** matrix, const char* string);

CMR_ERROR stringToCharMatrix(CMR* cmr, CMR_CHRMAT** matrix, const char* string);

#define ASSERT_CMR_CALL(x) \
  do \
  { \
    CMR_ERROR _error = (x); \
    ASSERT_EQ(_error, CMR_OKAY); \
  } \
  while (false)

#ifdef __cplusplus
}
#endif

#endif /* CMR_TEST_COMMON_H */

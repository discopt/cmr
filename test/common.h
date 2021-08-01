#ifndef TU_TEST_COMMON_H
#define TU_TEST_COMMON_H

#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

CMR_ERROR stringToDoubleMatrix(TU* tu, TU_DBLMAT** matrix, const char* string);

CMR_ERROR stringToIntMatrix(TU* tu, TU_INTMAT** matrix, const char* string);

CMR_ERROR stringToCharMatrix(TU* tu, TU_CHRMAT** matrix, const char* string);

#define ASSERT_TU_CALL(x) \
  ASSERT_FALSE(x)

#ifdef __cplusplus
}
#endif

#endif /* TU_TEST_COMMON_H */

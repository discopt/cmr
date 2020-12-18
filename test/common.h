#ifndef TU_TEST_COMMON_H
#define TU_TEST_COMMON_H

#include <tu/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

TU_ERROR stringToDoubleMatrix(TU* tu, TU_DBLMAT** matrix, const char* string);

TU_ERROR stringToIntMatrix(TU* tu, TU_INTMAT** matrix, const char* string);

TU_ERROR stringToCharMatrix(TU* tu, TU_CHRMAT** matrix, const char* string);

#ifdef __cplusplus
}
#endif

#endif /* TU_TEST_COMMON_H */

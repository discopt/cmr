#ifndef TU_TEST_COMMON_H
#define TU_TEST_COMMON_H

#include <tu/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

void stringToDoubleMatrix(TU* tu, TU_DBLMAT** matrix, const char* string);

void stringToIntMatrix(TU* tu, TU_INTMAT** matrix, const char* string);

void stringToCharMatrix(TU* tu, TU_CHRMAT** matrix, const char* string);

#ifdef __cplusplus
}
#endif

#endif /* TU_TEST_COMMON_H */

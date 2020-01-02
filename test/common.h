#ifndef TU_TEST_COMMON_H
#define TU_TEST_COMMON_H

#include <tu/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

TU_MATRIX_DOUBLE stringToMatrixDouble(const char* string);

TU_MATRIX_INT stringToMatrixInt(const char* string);

TU_MATRIX_CHAR stringToMatrixChar(const char* string);

#ifdef __cplusplus
}
#endif

#endif /* TU_TEST_COMMON_H */

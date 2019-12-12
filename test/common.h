#ifndef TU_TEST_COMMON_H
#define TU_TEST_COMMON_H

#include <tu/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

TU_SPARSE_DOUBLE stringToSparseDouble(const char* string);

TU_SPARSE_INT stringToSparseInt(const char* string);

TU_SPARSE_CHAR stringToSparseChar(const char* string);

#ifdef __cplusplus
}
#endif

#endif /* TU_TEST_COMMON_H */

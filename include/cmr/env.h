#ifndef CMR_ENV_H
#define CMR_ENV_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/config.h>
#include <cmr/export.h>
#include <stdbool.h>
#include <stdlib.h>

/**
 * \brief Type for return codes of library functions.
 **/

typedef enum
{
  CMR_OKAY = 0,         /**< No error. */
  CMR_ERROR_INPUT = 1,  /**< Bad user input. */
  CMR_ERROR_MEMORY = 2  /**< Error during (re)allocation. */
} CMR_ERROR;

/**
 * \brief Call wrapper for calls returning a \ref CMR_ERROR.
 */

#define TU_CALL(call) \
  do \
  { \
    CMR_ERROR _tu_error = call; \
    if (_tu_error) \
    { \
      if (_tu_error == CMR_ERROR_INPUT) \
        printf("User input error"); \
      else if (_tu_error == CMR_ERROR_MEMORY) \
        printf("Memory (re)allocation failed"); \
      else \
        printf("Unknown error"); \
      printf(" in %s:%d.\n", __FILE__, __LINE__); \
      return _tu_error; \
    } \
  } while (false)

struct TU_ENVIRONMENT;

/**
 * \brief Environment for computations
 *
 * Manages memory, threading, output and parameters.
 */

typedef struct TU_ENVIRONMENT TU;

/**
 * \brief Allocates and initializes a default \ref TU environment.
 *
 * It has default parameters and outputs to stdout.
 */

CMR_EXPORT
CMR_ERROR TUcreateEnvironment(
  TU** ptu /**< Pointer at which the \ref TU environment shall be allocated. */
);

/**
 * \brief Frees a \ref TU environment.
 */

CMR_EXPORT
CMR_ERROR TUfreeEnvironment(
  TU** ptu /**< Pointer to \ref TU environment. */
);

/**
 * \brief Allocates block memory for *\p ptr.
 *
 * Block memory shall be freed with \ref TUfreeBlock.
 * The size is determined automatically.
 */

#define TUallocBlock(tu, ptr) \
  _TUallocBlock(tu, (void**) ptr, sizeof(**ptr))

/**
 * \brief Carries out the allocation for \ref TUallocBlock.
 *
 * \note Use \ref TUallocBlock to allocate block memory.
 */

CMR_EXPORT
CMR_ERROR _TUallocBlock(TU* tu, void** ptr, size_t size);

/**
 * \brief Frees a block memory chunk allocated with \ref TUallocBlock.
 */

#define TUfreeBlock(tu, ptr) \
  _TUfreeBlock(tu, (void**) ptr, sizeof(**ptr))

/**
 * \brief Carries out the deallocation for \ref TUfreeBlock.
 *
 * \note Use \ref TUfreeBlock to free block memory.
 */

CMR_EXPORT
CMR_ERROR _TUfreeBlock(TU* tu, void** ptr, size_t size);

/**
 * \brief Allocates block memory for an array of chunks.
 *
 * The block memory shall be freed with \ref TUfreeBlockArray. Its size can be changed via
 * \ref TUreallocBlockArray.
 * The size of each chunk is determined automatically.
 */

#define TUallocBlockArray(tu, ptr, length) \
  _TUallocBlockArray(tu, (void**) ptr, sizeof(**ptr), length)

/**
 * \brief Carries out the allocation for \ref TUallocBlockArray.
 *
 * \note Use \ref TUallocBlockArray to allocate block memory.
 */

CMR_EXPORT
CMR_ERROR _TUallocBlockArray(TU* tu, void** ptr, size_t size, size_t length);

/**
 * \brief Reallocates block memory of an array of chunks.
 *
 * The block memory shall be freed with \ref TUfreeBlockArray.
 * The size of each chunk is determined automatically.
 */

#define TUreallocBlockArray(tu, ptr, length) \
  _TUreallocBlockArray(tu, (void**) ptr, sizeof(**ptr), length)

/**
 * \brief Carries out the reallocation for \ref TUreallocBlockArray.
 *
 * \note Use \ref TUreallocBlockArray to reallocate block memory.
 */

CMR_EXPORT
CMR_ERROR _TUreallocBlockArray(TU* tu, void** ptr, size_t size, size_t length);

/**
 * \brief Carries out the deallocation for \ref TUfreeBlockArray.
 *
 * \note Use \ref TUfreeBlockArray to free a block memory array.
 */

#define TUfreeBlockArray(tu, ptr) \
  _TUfreeBlockArray(tu, (void**) ptr)

CMR_EXPORT
CMR_ERROR _TUfreeBlockArray(TU* tu, void** ptr);

#ifdef __cplusplus
}
#endif

#endif /* CMR_ENV_H */


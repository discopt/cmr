#ifndef CMR_ENV_H
#define CMR_ENV_H

/**
 * \file env.h
 *
 * \author Matthias Walter
 *
 * \brief Basic functionality of the software library.
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/config.h>
#include <cmr/export.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* Macro for intended non-use of variables. */
#define CMR_UNUSED(x) (void)(x)


/**
 * \brief Type for return codes of library functions.
 **/

typedef enum
{
  CMR_OKAY = 0,               /**< No error. */
  CMR_ERROR_INPUT = 1,        /**< Bad user input. */
  CMR_ERROR_OUTPUT = 2,       /**< Error when writing user output. */
  CMR_ERROR_MEMORY = 3,       /**< Error during (re)allocation. */
  CMR_ERROR_INVALID = 4,      /**< Other invalid data. */
  CMR_ERROR_OVERFLOW = 5,     /**< Overflow in numerical computations. */
  CMR_ERROR_TIMEOUT = 6,      /**< Time limit exceeded. */
  CMR_ERROR_STRUCTURE = 7,    /**< Bad matrix structure. */
  CMR_ERROR_INCONSISTENT = 8, /**< Inconsistent pieces of input. */
} CMR_ERROR;

/**
 * \brief Call wrapper for calls returning a \ref CMR_ERROR.
 */

#define CMR_CALL(call) \
  do \
  { \
    CMR_ERROR _cmr_error = call; \
    if (_cmr_error) \
    { \
      if (_cmr_error == CMR_ERROR_INPUT) \
        fprintf(stderr, "User input error"); \
      else if (_cmr_error == CMR_ERROR_OUTPUT) \
        fprintf(stderr, "Error when writing user output"); \
      else if (_cmr_error == CMR_ERROR_MEMORY) \
        fprintf(stderr, "Memory (re)allocation failed"); \
      else if (_cmr_error == CMR_ERROR_INVALID) \
        fprintf(stderr, "Invalid input"); \
      else if (_cmr_error == CMR_ERROR_TIMEOUT) \
        fprintf(stderr, "Time limit exceeded"); \
      else if (_cmr_error == CMR_ERROR_OVERFLOW) \
        fprintf(stderr, "Integer overflow"); \
      else if (_cmr_error == CMR_ERROR_STRUCTURE) \
        fprintf(stderr, "Invalid matrix structure"); \
      else \
        fprintf(stderr, "Unknown error"); \
      fprintf(stderr, " in %s:%d.\n", __FILE__, __LINE__); \
      return _cmr_error; \
    } \
  } while (false)

struct CMR_ENVIRONMENT;

/**
 * \brief Environment for computations
 *
 * Manages memory, threading, output and parameters.
 */

typedef struct CMR_ENVIRONMENT CMR;

/**
 * \brief Allocates and initializes a default \ref CMR environment.
 *
 * It has default parameters and outputs to stdout.
 */

CMR_EXPORT
CMR_ERROR CMRcreateEnvironment(
  CMR** pcmr  /**< Pointer at which the \ref CMR environment shall be allocated. */
);

/**
 * \brief Frees a \ref CMR environment.
 */

CMR_EXPORT
CMR_ERROR CMRfreeEnvironment(
  CMR** pcmr  /**< Pointer to \ref CMR environment. */
);

/**
 * \brief Allocates block memory for *\p ptr.
 *
 * Block memory shall be freed with \ref CMRfreeBlock.
 * The size is determined automatically.
 */

#define CMRallocBlock(cmr, ptr) \
  _CMRallocBlock(cmr, (void**) ptr, sizeof(**ptr))

/**
 * \brief Carries out the allocation for \ref CMRallocBlock.
 *
 * \note Use \ref CMRallocBlock to allocate block memory.
 */

CMR_EXPORT
CMR_ERROR _CMRallocBlock(CMR* cmr, void** ptr, size_t size);

/**
 * \brief Frees a block memory chunk allocated with \ref CMRallocBlock.
 */

#define CMRfreeBlock(cmr, ptr) \
  _CMRfreeBlock(cmr, (void**) ptr, sizeof(**ptr))

/**
 * \brief Carries out the deallocation for \ref CMRfreeBlock.
 *
 * \note Use \ref CMRfreeBlock to free block memory.
 */

CMR_EXPORT
CMR_ERROR _CMRfreeBlock(CMR* cmr, void** ptr, size_t size);

/**
 * \brief Allocates block memory for an array of chunks.
 *
 * The block memory shall be freed with \ref CMRfreeBlockArray. Its size can be changed via
 * \ref CMRreallocBlockArray.
 * The size of each chunk is determined automatically.
 */

#define CMRallocBlockArray(cmr, ptr, length) \
  _CMRallocBlockArray(cmr, (void**) ptr, sizeof(**ptr), length)

/**
 * \brief Carries out the allocation for \ref CMRallocBlockArray.
 *
 * \note Use \ref CMRallocBlockArray to allocate block memory.
 */

CMR_EXPORT
CMR_ERROR _CMRallocBlockArray(CMR* cmr, void** ptr, size_t size, size_t length);

/**
 * \brief Reallocates block memory of an array of chunks.
 *
 * The block memory shall be freed with \ref CMRfreeBlockArray.
 * The size of each chunk is determined automatically.
 */

#define CMRreallocBlockArray(cmr, ptr, length) \
  _CMRreallocBlockArray(cmr, (void**) ptr, sizeof(**ptr), length)

/**
 * \brief Carries out the reallocation for \ref CMRreallocBlockArray.
 *
 * \note Use \ref CMRreallocBlockArray to reallocate block memory.
 */

CMR_EXPORT
CMR_ERROR _CMRreallocBlockArray(CMR* cmr, void** ptr, size_t size, size_t length);

/**
 * \brief Allocates block memory for an array of chunks and copies it from \p source.
 *
 * The block memory shall be freed with \ref CMRfreeBlockArray. Its size can be changed via
 * \ref CMRreallocBlockArray.
 * The size of each chunk is determined automatically.
 */

#define CMRduplicateBlockArray(cmr, ptr, length, source) \
  _CMRduplicateBlockArray(cmr, (void**) ptr, sizeof(**ptr), length, source)

/**
 * \brief Carries out the duplication for \ref CMRduplicateBlockArray.
 *
 * \note Use \ref CMRduplicateBlockArray to duplicate block memory.
 */

CMR_EXPORT
CMR_ERROR _CMRduplicateBlockArray(CMR* cmr, void** ptr, size_t size, size_t length, void* source);

/**
 * \brief Carries out the deallocation for \ref CMRfreeBlockArray.
 *
 * \note Use \ref CMRfreeBlockArray to free a block memory array.
 */

#define CMRfreeBlockArray(cmr, ptr) \
  _CMRfreeBlockArray(cmr, (void**) ptr)

CMR_EXPORT
CMR_ERROR _CMRfreeBlockArray(CMR* cmr, void** ptr);

/**
 * \brief Returns the error message (or \c NULL if no error with a message occured).
 */

CMR_EXPORT
char* CMRgetErrorMessage(
  CMR* cmr  /**< \ref CMR environment. */
);

/**
 * \brief Clears the error message.
 */

CMR_EXPORT
void CMRclearErrorMessage(
  CMR* cmr  /**< \ref CMR environment. */
);


#ifdef __cplusplus
}
#endif

#endif /* CMR_ENV_H */


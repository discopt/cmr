#include "io_internal.h"

#ifdef _WIN32

ssize_t getline(char** restrict lineptr, size_t* restrict n, FILE* restrict stream)
{
  assert(false);
}

#warning "_WIN32 is active"

#else

#warning "_WIN32 is NOT active"

#endif /* _WIN32 */

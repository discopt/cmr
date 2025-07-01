#include "io_internal.h"

#ifdef _WIN32

ssize_t getline(char** restrict lineptr, size_t* restrict n, FILE* restrict stream)
{
  assert(false);
}

#endif /* _WIN32 */

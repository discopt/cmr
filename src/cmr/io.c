#include "io_internal.h"

#include <errno.h>
#include <string.h>

#ifdef _WIN32

ssize_t getline(char** restrict lineptr, size_t* restrict n, FILE* restrict stream)
{
  assert(lineptr);
  assert(n);

  if (!lineptr || !n)
  {
    errno = EINVAL;
    return -1;
  }

  if (*lineptr == 0)
  {
    *n = 4;
    *lineptr = malloc((*n) * sizeof(char));
  }

  char* buffer = *lineptr;
  size_t remaining = *n;

  while (fgets(buffer, remaining, stream))
  {
    size_t read = strlen(buffer);

    if (read > 0 && buffer[read-1] == '\n')
      return strlen(*lineptr);
    else if (read > 0)
    {
      remaining = *n;
      (*n) *= 2;
      *lineptr = realloc(*lineptr, *n);
      // printf("Reallocated to %zu bytes: %p\n", *n, *lineptr);
      if (!*lineptr)
        return -1;
      buffer = &(*lineptr)[strlen(*lineptr)];
    }
  }

  return -1;
}

#endif /* _WIN32 */

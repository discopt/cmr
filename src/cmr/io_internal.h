#ifndef CMR_IO_INTERNAL_H
#define CMR_IO_INTERNAL_H

#include <stdio.h>

#include <cmr/env.h>


#ifdef _WIN32 /* getline is not available under Windows. */

ssize_t getline(char** restrict lineptr, size_t* restrict n, FILE* restrict stream);

#endif /* _WIN32 */


#endif /* CMR_IO_INTERNAL_H */

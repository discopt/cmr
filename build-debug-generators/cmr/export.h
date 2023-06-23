
#ifndef CMR_EXPORT_H
#define CMR_EXPORT_H

#ifdef CMR_STATIC_DEFINE
#  define CMR_EXPORT
#  define CMR_NO_EXPORT
#else
#  ifndef CMR_EXPORT
#    ifdef cmr_EXPORTS
        /* We are building this library */
#      define CMR_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define CMR_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef CMR_NO_EXPORT
#    define CMR_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef CMR_DEPRECATED
#  define CMR_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CMR_DEPRECATED_EXPORT
#  define CMR_DEPRECATED_EXPORT CMR_EXPORT CMR_DEPRECATED
#endif

#ifndef CMR_DEPRECATED_NO_EXPORT
#  define CMR_DEPRECATED_NO_EXPORT CMR_NO_EXPORT CMR_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CMR_NO_DEPRECATED
#    define CMR_NO_DEPRECATED
#  endif
#endif

#endif /* CMR_EXPORT_H */

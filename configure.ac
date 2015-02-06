AC_PREREQ(2.59)

AC_INIT([Unimodularity Library NoBoost], [1.2d], [https://github.com/xammy/unimodularity-test/wiki])

AM_INIT_AUTOMAKE([1.9 tar-pax])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADER([config.h])

AC_DEFINE([NDEBUG],[],[No-debug Mode])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.

# Checks for header files.

# Set language to C++, otherwise MAXOPT will check for CC compiler vendor and not CXX compiler vendor
AC_LANG(C++)

# Check for optimization and warnings
AX_COMPILER_VENDOR
AX_CXX_MAXOPT
AX_CXXFLAGS_WARN_ALL
AC_HEADER_ASSERT

AC_PROG_LIBTOOL

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


AC_CONFIG_FILES([
	Makefile
	src/Makefile
	])
AC_OUTPUT 

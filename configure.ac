AC_PREREQ(2.59)

AC_INIT([Total Unimodularity Library], [1.0], [xammy@xammy.homelinux.net])

AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADER([config.h])

#AC_DEFINE([NDEBUG],[],[No-debug Mode])

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

AC_PROG_LIBTOOL

AX_BOOST_BASE(1.33.1)
if test "x$want_boost" != "xyes"; then
	AC_MSG_ERROR([Sorry, but building without boost is not possible!])
fi


# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


AC_CONFIG_FILES([
	Makefile
	src/Makefile
	])
AC_OUTPUT 

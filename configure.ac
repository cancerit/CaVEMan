#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([CaVEMan], m4_esyscmd([tr -d '\n' < VERSION]), [cgpit@sanger.ac.uk], [$(PACKAGE_NAME)$(AC_PACKAGE_VERSION)], [https://github.com/cancerit/CaVEMan])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE(foreign subdir-objects)

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for libraries.
# Ensure zlib is available
AC_SEARCH_LIBS(gzbuffer, z, [], [AC_ERROR([zlib containing gzbuffer is required])])
AC_SEARCH_LIBS([tcase_create], [check_pic check, [], [AC_ERROR([check library not found])])
AC_SEARCH_LIBS(subunit_test_start, subunit)

PKG_CHECK_MODULES([zlib], [zlib > 1.2.3.3])

# Checks for header files.
AC_CHECK_HEADERS([fenv.h float.h inttypes.h limits.h stdint.h stdlib.h string.h unistd.h zlib.h math.h htslib/sam.h check.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_PROG_CC
AC_PROG_CXX

# Checks for pthread
AX_PTHREAD

AC_CANONICAL_HOST

case $host_os in
	*linux*)
		HOST_OS=linux
		;;
	*) ;;
esac
AM_CONDITIONAL([HOST_LINUX],[test x$HOST_OS = xlinux])

AC_FUNC_MALLOC
AC_CHECK_FUNCS([getcwd memset mkdir realpath rmdir strrchr])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

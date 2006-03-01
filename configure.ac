# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

AC_INIT([mincfit],[0.9],[Andrew Janke <a.janke@gmail.com>])
AC_CONFIG_SRCDIR([mincfit.c])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# check for volume_io
mni_REQUIRE_VOLUMEIO

# check for GSL
# AM_PATH_GSL(1.1)
AC_CHECK_LIB(m,main)
AC_CHECK_LIB(gslcblas,main)
AC_CHECK_LIB(gsl,main)

# Finish up by writing output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT(epm-header)

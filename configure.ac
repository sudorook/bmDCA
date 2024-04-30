dnl SPDX-FileCopyrightText: 2019 - 2022 sudorook <daemon@nullcodon.com>
dnl
dnl SPDX-License-Identifier: GPL-3.0-or-later
dnl
dnl This program is free software: you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by the Free
dnl Software Foundation, either version 3 of the License, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
dnl for more details.
dnl
dnl You should have received a copy of the GNU General Public License along
dnl with this program. If not, see <https://www.gnu.org/licenses/>.

AC_PREREQ([2.68])

AC_INIT([bmDCA],
        [0.9.5],
        [https://github.com/sudorook/bmDCA/issues],
        [bmDCA])

AC_CANONICAL_HOST

build_linux=no
build_macos=no
build_windows=no

# Detect the target system
case "${host_os}" in
  linux*)
    build_linux=yes
    ;;
  darwin*)
    build_macos=yes
    ;;
  cygwin*|mingw*|msys*)
    build_windows=yes
    ;;
  *)
    AC_MSG_ERROR(["'$host_os' is not supported"])
    ;;
esac

# Pass the conditionals to automake
AM_CONDITIONAL([LINUX], [test x$build_linux = xyes])
AM_CONDITIONAL([MACOS], [test x$build_macos = xyes])
AM_CONDITIONAL([WINDOWS], [test x$build_windows = xyes])

AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([enable])
AC_PREFIX_DEFAULT(/usr/local)
AM_SILENT_RULES([yes])

AC_LANG([C++])
AC_PROG_CXX

AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_cxx_compile_stdcxx.m4])
AX_CXX_COMPILE_STDCXX([17], [noext], [mandatory])

PKG_PROG_PKG_CONFIG

# Check for required libraries and headers.
PKG_CHECK_MODULES(
  [ARMADILLO],
  [armadillo],
  [],
  [AC_CHECK_HEADERS(
    [armadillo],
    [],
    [AC_MSG_ERROR("missing armadillo")],
    [])
  ])
AC_OPENMP
AC_CHECK_HEADERS(
  [stdio.h assert.h],
  [],
  [AC_MSG_ERROR("missing headers")],
  [])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_UINT16_T

AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT

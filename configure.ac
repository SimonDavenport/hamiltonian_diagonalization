##  This file generates the configure script for exact diagonalization programs

# Process this file with autoconf to produce a configure script.

AC_INIT(Exact Diagonalization, 1.0, simon.davenport2@gmail.com)
AC_CONFIG_SRCDIR(src/)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([subdir-objects])

# Option to use the mkl library
AC_ARG_ENABLE([intel_mkl],
  [AS_HELP_STRING([--enable-intel-mkl], [Set to yes to use intel mkl libray to speed up Lapack routines])])
  
# Option to use the blas library
AC_ARG_ENABLE([blas],
  [AS_HELP_STRING([--enable-blas], [Set to yes to add -lblas to library list])])

# Option to use the blas library
AC_ARG_ENABLE([lapack],
  [AS_HELP_STRING([--enable-lapack], [Set to yes to add -llapack to library list])])

# Option to use the math library
AC_ARG_ENABLE([math],
  [AS_HELP_STRING([--enable-math], [Set to yes to add -lm to library list])])

# Option to use a speed optimized map
AC_ARG_ENABLE([speed_optimization],
  [AS_HELP_STRING([--enable-speed-optimization], [Set to yes to use a speed optimized matrix  construction method (requires google dense hash)])])
  
# Option to use a memory optimized map
AC_ARG_ENABLE([memory_optimization],
  [AS_HELP_STRING([--enable-memory-optimization], [Set to yes to use a memory optimized matrix  construction method (requires google sparse hash)])])

# Extra arguments to configure google sparsehash
AC_ARG_WITH([google_hash_include],
            [AC_HELP_STRING([--with-google-hash-include=DIR],
                            [google sparsehash include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])

# Extra arguments to configure boost
AC_ARG_WITH([boost_include],
            [AC_HELP_STRING([--with-boost-include=DIR],
                            [boost include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
            
AC_ARG_WITH([boost_lib],
            [AC_HELP_STRING([--with-boost-lib=DIR],
                            [boost lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure fftw 
AC_ARG_WITH([fftw_include],
            [AC_HELP_STRING([--with-fftw-include=DIR],
                            [fftw include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
            
AC_ARG_WITH([fftw_lib],
            [AC_HELP_STRING([--with-fftw-lib=DIR],
                            [fftw lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure arpack
AC_ARG_WITH([arpack_lib],
            [AC_HELP_STRING([--with-arpack-lib=DIR],
                            [arpack lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure parpack
AC_ARG_WITH([parpack_lib],
            [AC_HELP_STRING([--with-parpack-lib=DIR],
                            [parpack lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure lapack
AC_ARG_WITH([lapack_lib],
            [AC_HELP_STRING([--with-lapack-lib=DIR],
                            [lapack lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure blas
AC_ARG_WITH([blas_lib],
            [AC_HELP_STRING([--with-blas-lib=DIR],
                            [blas lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure sqlite
AC_ARG_WITH([sqlite_include],
            [AC_HELP_STRING([--with-sqlite-include=DIR],
                            [sqlite3 include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
            
AC_ARG_WITH([sqlite_lib],
            [AC_HELP_STRING([--with-sqlite-lib=DIR],
                            [sqlite3 lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])
                                  
# Remove the default CXXFLAGS. 
CXXFLAGS="$CXXFLAGS"

# Add check for c++11 availability and set other compiler flags
CPPFLAGS="$CPPFLAGS -std=c++11 -O3 -Wall -Wfatal-errors"

# Check for a MPI c++ compiler, preferably one from intel
AC_PROG_CXX([mpicpc mpicxx mpicxxintel mpg++])

# Check for an intel fortran compiler, to be used for linking to intel mkl
# libraries only.
AS_IF([test "x$enable_intel_mkl" = "xyes"],[AC_PROG_FC([mpif90intel mpifort])])

# Otherwise set the linking compiler to be the c++ compiler
AS_IF([test "x$enable_intel_mkl" != "xyes"],[FC="$CXX"])

# Set the language to C++ in preparation for the further checks
AC_LANG(C++)

# Precompiler defs
AC_DEFINE(_DEBUG_,0)
AC_DEFINE(_ENABLE_ARPACK_,1)
AC_DEFINE(_USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_,1)
AC_DEFINE(_PYTHON_VERSION_,2)

AS_IF([test "x$enable_speed_optimization" = "xyes"],[AC_DEFINE(_SPEED_OPTIMIZED_MAP_,1) AC_DEFINE(_MEMORY_OPTIMIZED_MAP_,0) AC_CHECK_HEADERS(sparsehash/dense_hash_map)])

AS_IF([test "x$enable_memory_optimization" = "xyes"],[AC_DEFINE(_SPEED_OPTIMIZED_MAP_,0) AC_DEFINE(_MEMORY_OPTIMIZED_MAP_,1) AC_CHECK_HEADERS(sparsehash/sparse_hash_map)])

# Check for fftw3 header
AC_CHECK_HEADERS(fftw3.h)

# Check for sqlite3 header
AC_CHECK_HEADERS(sqlite3.h)

# Check for boost headers
AC_CHECK_HEADERS(boost/program_options.hpp)

# Add lm libraries to the LIB variabstdc++le
AS_IF([test "x$enable_math" = "xyes"],[LIBS="$LIBS -lm"])

# Add fortran compiler libraries to the LIB variable
AS_IF([test "x$enable_intel_mkl" = "xyes"],[LIBS="$LIBS -mkl -nofor_main -cxxlib"])

#   Add common library flags
LIBS="$LIBS -lparpack -larpack -lboost_program_options -lfftw3 -lsqlite3 -lgfortran"

# Add lapack libraries to the LIB variable
AS_IF([test "x$enable_lapack" = "xyes"],[LIBS="$LIBS -llapack"])

# Add blas libraries to the LIB variable
AS_IF([test "x$enable_blas" = "xyes"],[LIBS="$LIBS -lblas"])

# If not using a fortran compiler for linking, we need to ensure that
# the c++ compiler knows about fortran MPI functions when linking
AS_IF([test "x$enable_intel_mkl" != "xyes"],[LIBS="$LIBS -lmpi_f77 -lmpi_f90"])

# Check for ranlib, to generate the static library
AC_PROG_RANLIB

# Automatically generate Makefiles in each directory level
AC_OUTPUT(Makefile src/Makefile)

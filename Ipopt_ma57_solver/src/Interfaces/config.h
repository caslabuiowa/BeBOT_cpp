/* src/Common/config.h.  Generated from config.h.in by configure.  */
/* src/Common/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if the inexact linear solver option is included */
/* #undef BUILD_INEXACT */

/* Define to 1 if your C++ compiler doesn't accept -c and -o together. */
/* #undef CXX_NO_MINUS_C_MINUS_O */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define to 1 if your Fortran compiler doesn't accept -c and -o together. */
/* #undef F77_NO_MINUS_C_MINUS_O */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the <cfloat> header file. */
#define HAVE_CFLOAT 1

/* Define to 1 if you have the <cieeefp> header file. */
/* #undef HAVE_CIEEEFP */

/* Define to 1 if you have the <cmath> header file. */
#define HAVE_CMATH 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <float.h> header file. */
/* #undef HAVE_FLOAT_H */

/* Define to 1 if you have the <ieeefp.h> header file. */
/* #undef HAVE_IEEEFP_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <math.h> header file. */
/* #undef HAVE_MATH_H */

/* Define to 1 if you have the `MPI_Initialized' function. */
/* #undef HAVE_MPI_INITIALIZED */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vsnprintf' function. */
#define HAVE_VSNPRINTF 1

/* Define to 1 if windows.h is available. */
/* #undef HAVE_WINDOWS_H */

/* Define to 1 if you have the `_vsnprintf' function. */
/* #undef HAVE__VSNPRINTF */

/* Library Visibility Attribute */
/* #undef HSLLIB_EXPORT */

/* Library Visibility Attribute */
#define IPOPTAMPLINTERFACELIB_EXPORT 

/* Library Visibility Attribute */
#define IPOPTLIB_EXPORT 

/* Define to the debug sanity check level (0 is no test) */
#define IPOPT_CHECKLEVEL 0

/* Define to be the name of C-function for Inf check */
#define IPOPT_C_FINITE std::isfinite

/* Define to 1 if ASL is available. */
/* #undef IPOPT_HAS_ASL */

/* Define to 1 if function drand48 is available */
#define IPOPT_HAS_DRAND48 1

/* Define to 1 if function feenableexcept is available */
#define IPOPT_HAS_FEENABLEEXCEPT 1

/* Define to 1 if HSL is available. */
/* #undef IPOPT_HAS_HSL */

/* Define to 1 if the LAPACK package is available */
#define IPOPT_HAS_LAPACK 1

/* Define to 1 if the linear solver loader should be used to load libraries
   with linear solvers at runtime */
#define IPOPT_HAS_LINEARSOLVERLOADER 1

/* Define to 1 if Mumps is available. */
#define IPOPT_HAS_MUMPS 1

/* Define to 1 if Pardiso from MKL is available */
/* #undef IPOPT_HAS_PARDISO_MKL */

/* Define to 1 if function rand is available */
#define IPOPT_HAS_RAND 1

/* Define to 1 if the Spral package is available */
/* #undef IPOPT_HAS_SPRAL */

/* Define to 1 if function std::rand is available */
#define IPOPT_HAS_STD__RAND 1

/* Define to 1 if va_copy is available */
#define IPOPT_HAS_VA_COPY 1

/* Define to 1 if WSMP is available */
/* #undef IPOPT_HAS_WSMP */

/* Define to a macro mangling the given C identifier (in lower and upper
   case). */
/* #undef IPOPT_HSL_FUNC */

/* As IPOPT_HSL_FUNC, but for C identifiers containing underscores. */
/* #undef IPOPT_HSL_FUNC_ */

/* Define to 1 if Ipopt index type is int64_t */
/* #undef IPOPT_INT64 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case). */
#define IPOPT_LAPACK_FUNC(name,NAME) name ## _

/* As IPOPT_LAPACK_FUNC, but for C identifiers containing underscores. */
#define IPOPT_LAPACK_FUNC_(name,NAME) name ## _

/* Define to 1 if MPI should be initialized (finalized) when Ipopt library is
   (un)loaded. */
#define IPOPT_MPIINIT 1

/* Define to 1 if using single precision floating point */
/* #undef IPOPT_SINGLE */

/* Define to the debug verbosity level (0 is no output) */
#define IPOPT_VERBOSITY 0

/* Version number of project */
#define IPOPT_VERSION "3.14.13"

/* Major version number of project. */
#define IPOPT_VERSION_MAJOR 3

/* Minor version number of project. */
#define IPOPT_VERSION_MINOR 14

/* Release version number of project. */
#define IPOPT_VERSION_RELEASE 13

/* Define to a macro mangling the given C identifier (in lower and upper
   case). */
/* #undef IPOPT_WSMP_FUNC */

/* As IPOPT_WSMP_FUNC, but for C identifiers containing underscores. */
/* #undef IPOPT_WSMP_FUNC_ */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://github.com/coin-or/Ipopt/issues/new"

/* Define to the full name of this package. */
#define PACKAGE_NAME "Ipopt"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Ipopt 3.14.13"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "ipopt"

/* Define to the home page for this package. */
#define PACKAGE_URL "https://github.com/coin-or/Ipopt"

/* Define to the version of this package. */
#define PACKAGE_VERSION "3.14.13"

/* Name of Pardiso library from pardiso-project.org */
/* #undef PARDISO_LIB */

/* Library Visibility Attribute */
#define SIPOPTAMPLINTERFACELIB_EXPORT 

/* Library Visibility Attribute */
#define SIPOPTLIB_EXPORT 

/* The size of `int *', as computed by sizeof. */
#define SIZEOF_INT_P 8

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

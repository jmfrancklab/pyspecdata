/*
 * This file provides mappings for C functions callable by FORTRAN
 */

#ifdef ADD_UNDERSCORE
#define FORTRAN(n) n##_
#else
#define FORTRAN(n) n
#endif

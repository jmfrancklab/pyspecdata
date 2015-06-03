/* This is a kludge to work around the fact that VS FORTRAN won't let
   you write a line to standard output that does not end with a <CR>/<LF>
   unlike real programming languages.
*/
#include <stdio.h>
#include "fortrancall.h"

void FORTRAN(lprmpt)( void )
{
fputs( "NLSL> ", stdout );
fflush( stdout );
}

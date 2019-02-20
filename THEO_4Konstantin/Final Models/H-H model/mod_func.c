#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _fsin_reg();
extern void _fsquare_reg();
extern void _xtra_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," fsin.mod");
fprintf(stderr," fsquare.mod");
fprintf(stderr," xtra.mod");
fprintf(stderr, "\n");
    }
_fsin_reg();
_fsquare_reg();
_xtra_reg();
}

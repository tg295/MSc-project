#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _Nakpump_reg();
extern void _fsin_reg();
extern void _fsquare_reg();
extern void _hcn_reg();
extern void _ka_reg();
extern void _kdifl_reg();
extern void _kdr_reg();
extern void _kext_reg();
extern void _km_reg();
extern void _knatype_reg();
extern void _leak_reg();
extern void _nadifl_reg();
extern void _naext_reg();
extern void _nav1p7_reg();
extern void _nav1p8_reg();
extern void _nav1p9_reg();
extern void _xtra_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," Nakpump.mod");
fprintf(stderr," fsin.mod");
fprintf(stderr," fsquare.mod");
fprintf(stderr," hcn.mod");
fprintf(stderr," ka.mod");
fprintf(stderr," kdifl.mod");
fprintf(stderr," kdr.mod");
fprintf(stderr," kext.mod");
fprintf(stderr," km.mod");
fprintf(stderr," knatype.mod");
fprintf(stderr," leak.mod");
fprintf(stderr," nadifl.mod");
fprintf(stderr," naext.mod");
fprintf(stderr," nav1p7.mod");
fprintf(stderr," nav1p8.mod");
fprintf(stderr," nav1p9.mod");
fprintf(stderr," xtra.mod");
fprintf(stderr, "\n");
    }
_Nakpump_reg();
_fsin_reg();
_fsquare_reg();
_hcn_reg();
_ka_reg();
_kdifl_reg();
_kdr_reg();
_kext_reg();
_km_reg();
_knatype_reg();
_leak_reg();
_nadifl_reg();
_naext_reg();
_nav1p7_reg();
_nav1p8_reg();
_nav1p9_reg();
_xtra_reg();
}

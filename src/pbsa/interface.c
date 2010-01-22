#include "copyright_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include "prm.h"

// Author: Mengjuei Hsieh, University of California Irvine

typedef enum { FALSE, TRUE } boolean;

// epsin
// epsout
// istrng
// pbtemp
// dprob
// iprob
// accept
// maxitn
// fillratio
// space
// nbuffer
// nfocus
// fscale
// npbgrid
// arcres
// scalec
// cutres
// cutfd
// cutnb
// nsnbr
// nsnba
// phiout
// phiform
// npbverb
// use_rmin
// sprob
// vprob
// rhow_effect
// use_sav
// cavity_surften
// cavity_offset
// maxsph
// maxarc
// cutsa
// ndofd
// ndosas
// fmiccg
// ivalence
// laccept
// wsor
// lwsor
// pbkappa
// radinc
// expthresh
// offx
// offy
// offz
// sepbuf
// lmax
void main(){
	int    ipb;
	int    inp;
	int    smoothopt;
	int    radiopt;
	int    solvopt;
	int    npbopt;
	int    dbfopt;
	int    eneopt;
	int    frcopt;
	int    npopt;
	int    decompopt;
	int    bcopt;
	int    mpopt;
	double e_PBSA;

	get_e_pbsa(ipb,e_PBSA);
}

int get_e_pbsa(double e_PBSA){
	e_PBSA = 0;
	return 0;
}

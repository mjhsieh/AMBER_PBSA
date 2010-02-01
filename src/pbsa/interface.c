#include "copyright_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// Author: Mengjuei Hsieh, University of California Irvine

typedef	int	_INT_;
typedef	char	_STRING_;
typedef double	_REAL_;
extern FILE *nabout;

typedef enum { FALSE, TRUE } boolean;

typedef struct pbsa_opts {
	//Basic input options
	int    imin, ipb, inp;
	//Options to define the physical constants
	_REAL_ epsin, epsout, istrng, pbtemp, dprob, iprob, arcres;
	int    smoothopt, radiopt;
	//Options to select numerical procedures
	int    npbopt, solvopt, maxitn, nbuffer, fscale, npbgrid, nfocus;
	_REAL_ accept, fillratio, space, laccept, fmiccg, wsor, lwsor,
	       offx, offy, offz;
		//nfocus is undocumented
		//laccept is undocumented
		//fmiccg is undocumented
		//wsor is undocumented
		//lwsor is undocumented
		//offx is undocumented
		//offy is undocumented
		//offz is undocumented
	//Options to compute energy and forces
	int    bcopt, eneopt, frcopt, dbfopt, scalec, nsnbr, nsnba;
	_REAL_ cutres, cutfd, cutnb, cutsa;
		//cutsa is undocumented
	//Options for visualization and output
	int    phiout, phiform, npbverb;
	//Options to select a non-polar solvation treatment
	int    decompopt, use_rmin, use_sav, maxsph;
	_REAL_ sprob, vprob, rhow_effect, cavity_surften, cavity_offset;
	//Undocumented Options
	int    npopt, mpopt, ndofd, ndosas, lmax, maxarc;
	_REAL_ ivalence, radinc, expthresh, sepbuf;
} PBSA_OPTSSTRUCT_T;

typedef struct parm {
	char	ititl[81];
	_INT_ 	IfBox, Nmxrs, IfCap,
		 Natom,  Ntypes,  Nbonh,  Mbona,  Ntheth,  Mtheta, 
		 Nphih,  Mphia,  Nhparm, Nparm, Nnb, Nres,
		 Nbona,  Ntheta,  Nphia,  Numbnd,  Numang,  Nptra,
		 Natyp,  Nphb, Nat3, Ntype2d, Nttyp, Nspm, Iptres, Nspsol,
		 Ipatm, Natcap, Numextra;
	_STRING_ *AtomNames, *ResNames, *AtomSym, *AtomTree;
	_REAL_	*Charges, *Masses, *Rk, *Req, *Tk, *Teq, *Pk, *Pn, *Phase,
		 *Solty, *Cn1, *Cn2, *HB12, *HB10, *Rborn, *Fs;
	_REAL_	Box[4], Cutcap, Xcap, Ycap, Zcap, Fsmax;
	_INT_ 	*Iac, *Iblo, *Cno, *Ipres, *ExclAt, *TreeJoin, 
		 *AtomRes, *BondHAt1, *BondHAt2, *BondHNum, *BondAt1, *BondAt2, 
		 *BondNum, *AngleHAt1, *AngleHAt2, *AngleHAt3, *AngleHNum, 
		 *AngleAt1, *AngleAt2, *AngleAt3, *AngleNum, *DihHAt1, 
		 *DihHAt2, *DihHAt3, *DihHAt4, *DihHNum, *DihAt1, *DihAt2, 
		 *DihAt3, *DihAt4, *DihNum, *Boundary;
	_INT_	*N14pairs, *N14pairlist;
	_REAL_	*Gvdw;
} PARMSTRUCT_T;

PBSA_OPTSSTRUCT_T *pbsa_init(){
	PBSA_OPTSSTRUCT_T *myoptions;
	myoptions=(PBSA_OPTSSTRUCT_T *) malloc(sizeof(PBSA_OPTSSTRUCT_T));

	myoptions->imin			= 0;
	myoptions->ipb			= 0;
	myoptions->inp			= 0;
	
	myoptions->epsin		= 1.0;
	myoptions->epsout		= 80.0;
	myoptions->istrng		= 0.0;
	myoptions->pbtemp		= 300.0;
	myoptions->dprob		= 0.00;
	myoptions->iprob		= 2.00;
	myoptions->arcres		= 0.5;
	myoptions->smoothopt		= 0;
	myoptions->radiopt		= 1;

	myoptions->npbopt		= 0;
	myoptions->solvopt		= 2;
	myoptions->maxitn		= 100;
	myoptions->nbuffer		= 0;
	myoptions->fscale		= 8;
	myoptions->npbgrid		= 1;
	myoptions->nfocus		= 2;
	myoptions->accept		= 0.001;
	myoptions->fillratio		= 2.0;
	myoptions->space		= 0.5;
	myoptions->laccept		= 0.1;
	myoptions->fmiccg		= -0.30;
	myoptions->wsor			= 1.9;
	myoptions->lwsor		= 1.95;
	myoptions->offx			= 0.0;
	myoptions->offy			= 0.0;
	myoptions->offz			= 0.0;

	myoptions->bcopt		= 5;
	myoptions->eneopt		= -1;
	myoptions->frcopt		= 0;
	myoptions->dbfopt		= -1;
	myoptions->scalec		= 0;
	myoptions->nsnbr		= 1;
	myoptions->nsnba		= 1;
	myoptions->cutres		= 12.0;
	myoptions->cutfd		= 5.0;
	myoptions->cutnb		= 0.0;
	myoptions->cutsa		= 9.0;

	myoptions->phiout		= 0;
	myoptions->phiform		= 0;
	myoptions->npbverb		= 0;

	myoptions->decompopt		= 1;
	myoptions->use_rmin		= 0;
	myoptions->use_sav		= 0;
	myoptions->maxsph		= 400;

	myoptions->sprob		= 1.60;
	myoptions->vprob		= 1.28;
	myoptions->rhow_effect		= 1.0;
	myoptions->cavity_surften	= 0.04356;
	myoptions->cavity_offset	= 0.008;

	//undocumented
	myoptions->npopt		= 2;
	myoptions->mpopt		= 0;
	myoptions->ndofd		= 1;
	myoptions->ndosas		= 1;
	myoptions->lmax			= 80;
	myoptions->ivalence		= 1.0;
	myoptions->radinc		= myoptions->sprob*0.5;
	myoptions->maxarc		= 256;
	myoptions->expthresh		= 0.2;
	myoptions->sepbuf		= 4.0;

	return (myoptions); 
}


void main(){
	_REAL_ e_PBSA;

	get_e_pbsa(1,e_PBSA);
	printf("%.8f\n",e_PBSA);
}

int get_e_pbsa(int ipb, int inp, double e_PBSA){

	if ( ipb < 1 ) {
		printf("ipb should be either 1 or 2.\n");
		exit(1);
	}
        //initialize with default values
	PBSA_OPTSSTRUCT_T *pbsaopts;
	pbsaopts=pbsa_init();
	e_PBSA = 0;
	// not from human input

	free(pbsaopts);
	return 0;
}
/*
	int 	savbcopt[0] = 0;
	int 	savbcopt[1] = 0;
        int     dofd = 1;
	int     ntnbr = 1;
	int     ntnba = 1;
	boolean outphi = FALSE;
	boolean scalerf = FALSE;
	boolean srsas = TRUE;
	boolean pbverbose = FALSE;
	boolean pbprint = TRUE; 
	boolean pbgrid = TRUE;
	boolean pbinit = TRUE;
	boolean donpsa = TRUE;
 
   level = 1

   savxm(1) = 0
   savym(1) = 0
   savzm(1) = 0
   savxmym(1) = 0
   savxmymzm(1) = 0
   savh(1) = 0
   savxm(nfocus) = 0
   savym(nfocus) = 0
   savzm(nfocus) = 0
   savxmym(nfocus) = 0
   savxmymzm(nfocus) = 0
   savh(nfocus) = 0
  
   nsaslag = 100
   lastp = 0
   pbgamma_int = 1.0
   pbgamma_ext = 65.0
	myoptions->maxarcdot	= 1500;
*/

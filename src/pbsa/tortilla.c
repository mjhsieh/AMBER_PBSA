#include "copyright_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// Author: Mengjuei Hsieh, University of California Irvine

typedef	int	INT_T;
typedef	char	_STRING_;
typedef double	REAL_T;

typedef enum { FALSE, TRUE } boolean;

typedef struct pbsa_opts {
	//Basic input options
	INT_T  imin, ipb, inp;
	//Options to define the physical constants
	REAL_T epsin, epsout, istrng, pbtemp, dprob, iprob, arcres;
	INT_T  smoothopt, radiopt;
	//Options to select numerical procedures
	INT_T  npbopt, solvopt, maxitn, nbuffer, fscale, npbgrid, nfocus;
	REAL_T accept, fillratio, space, laccept, fmiccg, wsor, lwsor,
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
	INT_T  bcopt, eneopt, frcopt, dbfopt, scalec, nsnbr, nsnba;
	REAL_T cutres, cutfd, cutnb, cutsa;
		//cutsa is undocumented
	//Options for visualization and output
	INT_T  phiout, phiform, npbverb;
	//Options to select a non-polar solvation treatment
	INT_T  decompopt, use_rmin, use_sav, maxsph;
	REAL_T sprob, vprob, rhow_effect, cavity_surften, cavity_offset;
	//Undocumented Options
	INT_T  npopt, mpopt, ndofd, ndosas, lmax, maxarc;
	REAL_T ivalence, radinc, expthresh, sepbuf;
} PBSA_OPTSSTRUCT_T;

//extern FILE *nabout;

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


REAL_T epbsa(int ipb, REAL_T fillratio, INT_T natom,
	       int nres, int ntype, int *ipres, int *iac,
	       int ico, REAL_T *x){
        REAL_T e_PBSA;
	//checking passed options
	if ( ipb != 1 && ipb !=2 ) {
		printf("ipb should be either 1 or 2.\n");
		exit(1);
	/*
	} else if ( nfocus != 1 && nfocus != 2) {
		printf("nfocus should be either 1 or 2.\n");
		exit(1);
	*/
	} else if ( fillratio < 1 ) {
		printf("fillratio should be at least 1.\n");
		exit(1);
	}
        //initialize with default values
	PBSA_OPTSSTRUCT_T *pbsaopts;
	pbsaopts=pbsa_init();
	e_PBSA = 0;
	//overridding defaults
	pbsaopts->ipb=ipb;
	pbsaopts->fillratio=fillratio;
	// not from human input
	// int ix(i02) RESIDUE_POINTER
	// int ix(i04) ATOM_TYPE_INDEX
	// int ix(i06) INDEX TO N-B TYPE
        //pb_force_(natom,nres,ntypes,ix(i02),ix(i04),ix(i06),ix(i10),cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
	//np_force_(natom,nres,ntypes,ix(i02),ix(i04),ix(i06),cn1,cn2,x,f,esurf,edisp)
	free(pbsaopts);
	return (e_PBSA);
}
/*
	int 	savbcopt[0] = 0;
	int 	savbcopt[1] = 0;
        INT_T   dofd = 1;
	INT_T   ntnbr = 1;
	INT_T   ntnba = 1;
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
   pbgammaINT_T = 1.0
   pbgamma_ext = 65.0
	myoptions->maxarcdot	= 1500;
*/

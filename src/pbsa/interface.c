#include "copyright_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// Author: Mengjuei Hsieh, University of California Irvine

typedef	int	_INT_;
typedef	char	_STRING_;
typedef double	_REAL_;

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
PARMSTRUCT_T *rdparm();
PARMSTRUCT_T *myparm;

extern FILE *nabout;

void main(){
	_REAL_ e_PBSA;
	_INT_  ipb=1;
	_REAL_ fillratio=4;


	get_e_pbsa(ipb,fillratio,e_PBSA);
	printf("%.8f\n",e_PBSA);
}

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


int get_e_pbsa(int ipb, _REAL_ fillratio,
		_REAL_ e_PBSA){

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

	free(pbsaopts);
	return (0);
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



int trim_rspc(char *mystring){
    int i,len;
    boolean bool_end=FALSE;
    len=strlen(mystring);
    i=len-1;
    while (i>=-1 && !bool_end){
        if(mystring[i]==' '){
            mystring[i]='\0';
        }else{
            bool_end=TRUE;
        }
        i--;
    }
}

int open_parm_(char *prm_path){
    char prm_fn[80];

    strcpy(prm_fn, prm_path);
    trim_rspc(prm_fn);

    if ((nabout = fopen("/dev/null","w"))==NULL){
        printf("having problem opening nabout\n");
        exit(1);
    }

    myparm = rdparm(prm_fn);

    fclose(nabout);
    return (0);
}

int rdprm1_(char *title,int *natom,int *ntypes,int *nbonh,int *mbona,
		int *ntheth,int *mtheta,int *nphih,int *mphia,
		int *nhparm,int *nparm,int *nnb,int *nres,int *nbona,
		int *ntheta,int *nphia,int *numbnd,int *numang,
		int *nptra,int *natyp,int *nphb,int *ifbox,int *nmxrs,
		int *ifcap,int *numextra,int *ncopy){
	int i;
	//strcpy(title,myparm->ititl);
	for (i=0;i<80;i++){
		title[i]=myparm->ititl[i];
	}
	*natom=myparm->Natom;
	*ntypes=myparm->Ntypes;
	*nbonh=myparm->Nbonh;
	*mbona=myparm->Mbona;
	*ntheth=myparm->Ntheth;
	*mtheta=myparm->Mtheta;
	*nphih=myparm->Nphih;
	*mphia=myparm->Mphia;
	*nhparm=myparm->Nhparm;
	*nparm=myparm->Nparm;
	*nnb=myparm->Nnb;
	*nres=myparm->Nres;
	*nbona=myparm->Nbona;
	*ntheta=myparm->Ntheta;
	*nphia=myparm->Nphia;
	*numbnd=myparm->Numbnd;
	*numang=myparm->Numang;
	*nptra=myparm->Nptra;
	*natyp=myparm->Natyp;
	*nphb=myparm->Nphb;
	*ifbox=myparm->IfBox;
	*nmxrs=myparm->Nmxrs;
	*ifcap=myparm->IfCap;
	*numextra=myparm->Numextra;
}

int rdprm2_(char *myAtomName,double *myCharge,double *myMass,int *myAtmType,
    int *myNXcldAtm,int *ntype,int *myNBPrm,char *myResNam,int *myResPtr,
    double *myRk,double *myReq,double *myTk,double *myTeq,double *myPk,
    double *myPn,double *myPhase,double *mySolty,
    int *nttyp,double *myCn1,double *myCn2,
    int *myBondHAt1,int *myBondHAt2,int *myBondHNum,
    int *myBondAt1,int *myBondAt2,int *myBondNum,
    int *myAngleHAt1,int *myAngleHAt2,int *myAngleHAt3,int *myAngleHNum,
    int *myAngleAt1,int *myAngleAt2,int *myAngleAt3,int *myAngleNum,
    int *myDihHAt1,int *myDihHAt2,int *myDihHAt3,int *myDihHAt4,int *myDihHNum,
    int *myDihAt1,int *myDihAt2,int *myDihAt3,int *myDihAt4,int *myDihNum,
    int *myExclAt,int *myHB12,int *myHB10,
    char *myAtomSym,char *myAtomTree,int *myTreeJoin,
    int *myAtomRes){
        int i;
	int natom=myparm->Natom;
	for (i=0;i<natom*4;i++){
		myAtomName[i]=myparm->AtomNames[i];
		myAtomSym[i] =myparm->AtomSym[i];
		myAtomTree[i]=myparm->AtomTree[i];
	}
	for (i=0;i<natom;i++){
		myCharge[i]  =myparm->Charges[i];
		myMass[i]    =myparm->Masses[i];
		myAtmType[i] =myparm->Iac[i];
		myNXcldAtm[i]=myparm->Iblo[i];
		myTreeJoin[i]=myparm->TreeJoin[i];
		myAtomRes[i] =myparm->AtomRes[i];
	}
	for (i=0;i<*ntype;i++){
		myNBPrm[i]=myparm->Cno[i];
	}
	int nres=myparm->Nres;
	for (i=0;i<nres*4;i++){
		myResNam[i]=myparm->ResNames[i];
	}
	for (i=0;i<nres;i++){
		myResPtr[i]=myparm->Ipres[i];
	}
	int numbnd=myparm->Numbnd;
	for (i=0;i<numbnd;i++){
		myRk[i] =myparm->Rk[i];
		myReq[i]=myparm->Req[i];
	}
	int numang=myparm->Numang;
	for (i=0;i<numang;i++){
		myTk[i] =myparm->Tk[i];
		myTeq[i]=myparm->Teq[i];
	}
	int nptra=myparm->Nptra;
	for (i=0;i<nptra;i++){
		myPk[i]=myparm->Pk[i];
		myPn[i]=myparm->Pn[i];
		myPhase[i]=myparm->Phase[i];
	}
	int natyp=myparm->Natyp;
	for (i=0;i<natyp;i++){
		mySolty[i]=myparm->Solty[i];
	}
	for (i=0;i<*nttyp;i++){
		myCn1[i]=myparm->Cn1[i];
		myCn2[i]=myparm->Cn2[i];
	}
	int nbonh=myparm->Nbonh;
	for (i=0;i<nbonh;i++){
		myBondHAt1[i]=myparm->BondHAt1[i];
		myBondHAt2[i]=myparm->BondHAt2[i];
		myBondHNum[i]=myparm->BondHNum[i];
	}
	int nbona=myparm->Nbona;
	for (i=0;i<nbona;i++){
		myBondAt1[i]=myparm->BondAt1[i];
		myBondAt2[i]=myparm->BondAt2[i];
		myBondNum[i]=myparm->BondNum[i];
	}
	int ntheth=myparm->Ntheth;
	for (i=0;i<ntheth;i++){
		myAngleHAt1[i]=myparm->AngleHAt1[i];
		myAngleHAt2[i]=myparm->AngleHAt2[i];
		myAngleHAt3[i]=myparm->AngleHAt3[i];
		myAngleHNum[i]=myparm->AngleHNum[i];
	}
	int ntheta=myparm->Ntheta;
	for (i=0;i<ntheta;i++){
		myAngleAt1[i]=myparm->AngleAt1[i];
		myAngleAt2[i]=myparm->AngleAt2[i];
		myAngleAt3[i]=myparm->AngleAt3[i];
		myAngleNum[i]=myparm->AngleNum[i];
	}
	int nphih=myparm->Nphih;
	for (i=0;i<nphih;i++){
		myDihHAt1[i]=myparm->DihHAt1[i];
		myDihHAt2[i]=myparm->DihHAt2[i];
		myDihHAt3[i]=myparm->DihHAt3[i];
		myDihHAt4[i]=myparm->DihHAt4[i];
		myDihHNum[i]=myparm->DihHNum[i];
	}
	int nphia=myparm->Nphia;
	for (i=0;i<nphia;i++){
		myDihAt1[i]=myparm->DihAt1[i];
		myDihAt2[i]=myparm->DihAt2[i];
		myDihAt3[i]=myparm->DihAt3[i];
		myDihAt4[i]=myparm->DihAt4[i];
		myDihNum[i]=myparm->DihNum[i];
	}
	int nnb=myparm->Nnb;
	for (i=0;i<nnb;i++){
		myExclAt[i]=myparm->ExclAt[i];
	}
	int nphb=myparm->Nphb;
	for (i=0;i<nphb;i++){
		myHB12[i]=myparm->HB12[i];
		myHB10[i]=myparm->HB10[i];
		//myHBCUT[i]=myparm->H[i];//not in this struct
	}
}

int rdprm3_(int *myIptres,int *myNspm,int *myNspsol,int *myIpatm){
	*myIptres=myparm->Iptres;
	*myNspm=myparm->Nspm;
	*myNspsol=myparm->Nspsol;
	*myIpatm=myparm->Ipatm;
}

int rdprm4_(int *myNatcap,double *myCutcap,double *myXcap,double *myYcap,
		double *myZcap){
	*myNatcap=myparm->Natcap;
	*myCutcap=myparm->Cutcap;
	*myXcap=myparm->Xcap;
	*myYcap=myparm->Ycap;
	*myZcap=myparm->Zcap;
}

int rdprm5_(double *myRborn,double *myFs){
        int i;
	int natom=myparm->Natom;
	for (i=0;i<natom;i++){
		myRborn[i]=myparm->Rborn[i];
		myFs[i]=myparm->Fs[i];
	}
}

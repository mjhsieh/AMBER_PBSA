#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "prm.h"

// Author: Mengjuei Hsieh, University of California Irvine

#include "is_copyrightc.h"

typedef enum { FALSE, TRUE } boolean;

PARMSTRUCT_T *rdparm();
PARMSTRUCT_T *myparm;

FILE *nabout;

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

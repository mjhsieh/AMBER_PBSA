#include "copyright_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "prm.h"

// Author: Mengjuei Hsieh, University of California Irvine

typedef enum { FALSE, TRUE } boolean;

//extern FILE *nabout;
void pb_read_(int*,REAL_T*);
void pb_init_(int*,INT_T*,INT_T*,INT_T*,INT_T*,INT_T*,
	     INT_T*,INT_T*,INT_T*,INT_T*,INT_T*,
	     INT_T*,INT_T*,INT_T*,INT_T*,
	     STRING_T*,STRING_T*,STRING_T*,REAL_T*,
	     REAL_T*);
void mypb_force_(INT_T*,INT_T*,INT_T*,
		INT_T*,INT_T*,INT_T*,INT_T*,
		REAL_T*,REAL_T*,REAL_T*,REAL_T*,REAL_T*,REAL_T*,
		REAL_T*,REAL_T*,REAL_T*,REAL_T*);
void pb_free_();

// int ix(i02) RESIDUE_POINTER
// int ix(i04) ATOM_TYPE_INDEX
// int ix(i06) INDEX TO N-B TYPE
// int ix(i10) EXCLUDED ATOM LIST
// pb_force_(natom,       nres,       ntypes,       ix(i02)
//           ix(i04),     ix(i06),    ix(i10),      cn1,
//           cn2,         xx(l15),    x,            f,
//           evdw,        eelt,       epol                             )
//REAL_T epbsa(int ipb, REAL_T fillratio,
//	     INT_T natom, INT_T nres, INT_T ntypes, INT_T *ipres,
//	     INT_T *iac,  INT_T *ico, INT_T *exclat,REAL_T *cn1,
//	     REAL_T *cn2, REAL_T *cg, REAL_T *x,    REAL_T *f){
REAL_T epbsa(int ipb, REAL_T fillratio, PARMSTRUCT_T *prm, REAL_T *x,
	     REAL_T *grad, REAL_T *evdw, REAL_T *eelt, REAL_T *esurf,
	     REAL_T *edisp){
	int i, ifcap;
	ifcap = 0;
        REAL_T e_PBSA=0;
	REAL_T *f;
	//checking passed options
	if ( ipb != 1 && ipb !=2 ) {
		printf("ipb should be either 1 or 2.\n"); exit(1);
	} else if ( fillratio < 1 ) {
		printf("fillratio should be at least 1.\n"); exit(1);
	}

	pb_read_(&ipb,&fillratio);
//call pb_init(ifcap,natom,nres,ntypes,nbonh,nbona,
//             ix(i02),ix(i04),ix(i06),ix(i08),ix(i10),
//             ix(iibh),ix(ijbh),ix(iiba),ix(ijba),
//             ih(m02),ih(m04),ih(m06),x(l15),
//             x(l97))
        //Arrays: Ipres,Iac,Cno,ExclAt,Cn1,Cn2,x,grad
	prepb_init_(prm->Cn1,prm->Cn2,&prm->Nttyp);
        pb_init_(&ifcap,&prm->Natom,&prm->Nres,&prm->Ntypes,&prm->Nbonh,&prm->Nbona,
	     prm->Ipres,prm->Iac,prm->Cno,prm->Iblo,prm->ExclAt,
	     prm->BondHAt1,prm->BondHAt2,prm->BondAt1,prm->BondAt2,
	     prm->ResNames,prm->AtomNames,prm->AtomSym,prm->Charges,
	     prm->Rborn);

	f = (REAL_T *) malloc(sizeof(REAL_T)*3*prm->Natom);

	mypb_force_(&prm->Natom,&prm->Nres,&prm->Ntypes,
	     prm->Ipres,prm->Iac,prm->Cno,prm->ExclAt,
	     prm->Cn1,prm->Cn2,prm->Charges,x,f,&e_PBSA,
	     evdw,eelt,esurf,edisp);

	for (i=0;i<3*prm->Natom;i++){
		grad[i]=grad[i]-f[i];
	}

	free(f);
	pb_free_();
	return (e_PBSA);
}

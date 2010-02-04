#include "copyright_c.h"
// Description:
//	A program to test C interface for libpbsa
// Author:
// 	Mengjuei Hsieh, University of California Irvine
//
#include <stdio.h>
#include <stdlib.h>
#include "prm.h"

FILE *nabout;
PARMSTRUCT_T *prm;

void getx_(int*,int*,REAL_T*);
REAL_T epbsa(int,REAL_T,INT_T,INT_T,INT_T,
	     INT_T*,INT_T*,INT_T*,INT_T*,REAL_T*,
	     REAL_T*,REAL_T*,REAL_T*,REAL_T*);

int main(){
    int    pbsa;
    double fillratio;
    REAL_T e_pb;
    char   *prm_fpath="unittest.prmtop";
    int    filenum,i;
    REAL_T *x,*grad;

    pbsa = 1;
    fillratio = 4.0;
    e_pb = 0.0;
    filenum = 9;

    if ((nabout = fopen("/dev/null","w"))==NULL){
        printf("having problem opening nabout\n");
        exit (1);
    }

    prm = (PARMSTRUCT_T *) malloc(sizeof(PARMSTRUCT_T));
    prm = rdparm(prm_fpath); fclose(nabout);

    x = (REAL_T *) malloc(sizeof(REAL_T)*3*prm->Natom);
    grad = (REAL_T *) malloc(sizeof(REAL_T)*3*prm->Natom);
    getx_(&prm->Natom,&filenum,x); //to import x from inpcrd

    //Arrays: Ipres,Iac,Cno,ExclAt,Cn1,Cn2,x,grad
    e_pb = epbsa(pbsa,fillratio,prm->Natom,prm->Nres,prm->Ntypes,
                prm->Ipres,prm->Iac,prm->Cno,prm->ExclAt,prm->Cn1,
                prm->Cn2,prm->Charges,x,grad);
    printf("EPB = %.4f\n",e_pb);

    free(x); free(grad); free(prm);
    return (0);
}

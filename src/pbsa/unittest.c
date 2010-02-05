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

void private_getx_(int*,int*,REAL_T*);
REAL_T epbsa(int,REAL_T,PARMSTRUCT_T*,REAL_T*,REAL_T*,
	    REAL_T*,REAL_T*,REAL_T*,REAL_T*);

int main(){
    int    pbsa;
    double fillratio;
    REAL_T e_pb,evdw,eelt,esurf,edisp;
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
    private_getx_(&prm->Natom,&filenum,x); //to import x from inpcrd

    e_pb = epbsa(pbsa,fillratio,prm,x,grad,
		 &evdw,&eelt,&esurf,&edisp);
    printf("VDWAALS = %12.4f\n",evdw);
    printf("EPB     = %12.4f\n",e_pb);
    printf("EELEC   = %12.4f\n",eelt);
    printf("ECAVITY = %12.4f\n",esurf);
    printf("EDISPER = %12.4f\n",edisp);

    free(x); free(grad); free(prm);
    return (0);
}

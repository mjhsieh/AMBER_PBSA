typedef struct pbopts {
    REAL_T epsin, epsout, istrng, pbtemp, dprob, iprob, accept,
           fillratio, space, arcres, cutres, cutfd,
           cutnb, sprob, vprob, rhow_effect, cavity_surften,
           cavity_offset, cutsa, fmiccg, ivalence, laccept, wsor,
           lwsor, pbkappa, radinc, expthresh, offx, offy, offz,
           sepbuf;
    int    smoothopt, radiopt, npbopt, solvopt, maxitn, nbuffer,
           nfocus, fscale, npbgrid, dbfopt, bcopt, scalec, eneopt,
           frcopt, nsnbr, phiout, phiform, npbverb, npopt,
           decompopt, use_rmin, use_sav, maxsph, maxarc, ndofd,
           ndosas, mpopt, lmax, inp;
} PBOPTSTRUCT_T;

PBOPTSTRUCT_T* pboptinit();

REAL_T epbsa(int ipb, PBOPTSTRUCT_T *opt, PARMSTRUCT_T *prm, REAL_T *x,
	REAL_T *grad, REAL_T *evdw, REAL_T *eelt, REAL_T *esurf,
	REAL_T *edisp);


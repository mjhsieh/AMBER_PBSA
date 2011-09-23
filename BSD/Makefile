#
#************************************************************************
# 
#************************************************************************
#
include ../config.h

CONFIG_FILE=../config.h
CONFIG_COMMAND=configure
PARALLEL_TOKEN=DMPI

# PBSA Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

SRC=\
	dsvdc.f  pbsa.f pb_read.f pb_write.f getcoor.f \
	runmd.f runmin.f force.f \
	pb_init.f sa_driver.f np_force.f pb_force.f pb_fddrv.f pb_iimdrv.f pb_nlsolver.f \
	pb_exmol.f pb_mpfrc.f pb_direct.f pb_list.f timer.f \
	egb.f ene.f locmem.f myopen.f rdparm.f pb_lsolver.f \
	decomp.f rgroup.f rfree.f debug.f \
	coed20.f  indexg.f  irre32.f  phidev.f  project.f  qld.f \
        coed6.f   irre31.f  miniop.f  prodis.f  qint.f     wint.f \
        transf.f curv.f GrToPr.f jumps.f problem.f \
        IIM.f regular.f ctime.f gmres.f \
	daxpy.f  ddot.f    dheqr.f  dpigmr.f  dscal.f   dslui2.f  dxlcal.f \
	dchkw.f  dgmres.f  dnrm2.f  drlcal.f  dsilus.f  dslui.f   isdgmr.f \
	dcopy.f  dhels.f   dorth.f  ds2y.f    dslugm.f  dsmv.f    qs2i1d.f \
	d1mach.f  fdump.f   i1mach.f  j4save.f  xercnt.f  xerhlt.f  xermsg.f  xerprn.f  xersve.f  xgetua.f \
        bicg.f dbcg.f dslui4.f dsluti.f dsmtv.f dslubc.f isdbcg.f pb_chkfrc.f svdcmp.f svbksb.f pythag.f

OBJ=\
	dsvdc.o  pbsa.o pb_read.o pb_write.o getcoor.o \
	runmd.o runmin.o force.o \
	pb_init.o sa_driver.o np_force.o pb_force.o pb_fddrv.o pb_iimdrv.o pb_nlsolver.o \
	pb_exmol.o pb_mpfrc.o pb_direct.o pb_list.o timer.o \
	egb.o ene.o locmem.o myopen.o rdparm.o pb_lsolver.o \
	decomp.o rgroup.o rfree.o debug.o \
	coed20.o  indexg.o  irre32.o  phidev.o  project.o  qld.o \
	coed6.o   irre31.o  miniop.o  prodis.o  qint.o     wint.o \
	transf.o curv.o GrToPr.o jumps.o problem.o \
	IIM.o regular.o ctime.o gmres.o \
	daxpy.o  ddot.o    dheqr.o  dpigmr.o  dscal.o   dslui2.o  dxlcal.o \
	dchkw.o  dgmres.o  dnrm2.o  drlcal.o  dsilus.o  dslui.o   isdgmr.o \
	dcopy.o  dhels.o   dorth.o  ds2y.o    dslugm.o  dsmv.o    qs2i1d.o \
	d1mach.o  fdump.o   i1mach.o  j4save.o  xercnt.o  xerhlt.o  xermsg.o  \
        xerprn.o  xersve.o  xgetua.o \
	bicg.o dbcg.o dslui4.o dsluti.o dsmtv.o dslubc.o isdbcg.o pb_chkfrc.o svdcmp.o svbksb.o pythag.o

LIBPBSAOBJS= \
	interface.o pbsaapi.o   myopen.o    pb_read.LIBPBSA.o pb_force.o \
	pb_fddrv.o  pb_iimdrv.o sa_driver.o pb_nlsolver.o     pb_lsolver.o \
	dsvdc.o     timer.o     pb_init.o   pb_exmol.o        np_force.o \
	pb_direct.o pb_list.o   pb_write.o  pb_mpfrc.o        decomp.o \
	rgroup.o    rfree.o     coed20.o    indexg.o          irre32.o \
	phidev.o    project.o   qld.o       coed6.o           irre31.o \
	miniop.o    prodis.o    qint.o      wint.o            transf.o \
	curv.o      GrToPr.o    jumps.o     problem.o         IIM.o \
	regular.o         ctime.o \
        gmres.o     bicg.o      project.o   transf.o          dslubc.o \
        dslugm.o    curv.o      irre31.o    irre32.o          \
	problem.o   wint.o      dchkw.o     ds2y.o            GrToPr.o \
	dgmres.o    dslui.o     dsluti.o \
	dsmv.o      dsmtv.o     dsilus.o    jumps.o           dbcg.o \
	qld.o       phidev.o    dpigmr.o    dslui2.o          dslui4.o \
	isdbcg.o    qs2i1d.o    xermsg.o    d1mach.o          xercnt.o \
	xerhlt.o    drlcal.o    j4save.o    isdgmr.o          fdump.o \
	xerprn.o    dorth.o     dhels.o     xersve.o          dheqr.o \
	i1mach.o    xgetua.o    dxlcal.o    \
        svdcmp.o    svbksb.o    pythag.o

SANDERPBSAOBJS= dsvdc.o sa_driver.SANDER.o pb_nlsolver.o \
        pb_lsolver.o pb_force.SANDER.o np_force.SANDER.o \
        pb_exmol.SANDER.o pb_mpfrc.SANDER.o pb_direct.SANDER.o \
        pb_list.SANDER.o pb_read.SANDER.o pb_write.SANDER.o pb_init.SANDER.o \
	pb_fddrv.o pb_iimdrv.o timer.SANDER.o \
	coed20.o  indexg.o  irre32.o  phidev.o  project.o  qld.o \
	coed6.o   irre31.o  miniop.o  prodis.o  qint.o     wint.o \
	transf.o curv.o GrToPr.o jumps.o problem.o \
	IIM.o    regular.o ctime.o \
	gmres.o  daxpy.o  ddot.o   dheqr.o  dpigmr.o dscal.o  dslui2.o \
        dxlcal.o dchkw.o  dgmres.o dnrm2.o  drlcal.o dsilus.o dslui.o  \
        isdgmr.o dcopy.o  dhels.o  dorth.o  ds2y.o   dslugm.o dsmv.o   \
        qs2i1d.o d1mach.o fdump.o  i1mach.o j4save.o xercnt.o xerhlt.o \
        xermsg.o xerprn.o xersve.o xgetua.o \
        bicg.o   dbcg.o   dslui4.o dsluti.o dsmtv.o  dslubc.o isdbcg.o \
	svdcmp.o svbksb.o pythag.o

install: install.$(INSTALLTYPE)

install.serial: serial
	/bin/mv pbsa$(SFX) $(BINDIR)

install.parallel: parallel
	/bin/mv pbsa.MPI$(SFX) $(BINDIR)

serial: configured_serial pbsa$(SFX)

parallel: configured_parallel pbsa.MPI$(SFX)

uninstall: uninstall.$(INSTALLTYPE)

uninstall.serial: configured_serial
	-rm -f $(BINDIR)/pbsa$(SFX)
	#-rm -f $(BINDIR)/simplepbsa$(SFX)

uninstall.parallel: configured_parallel
	-rm -f $(BINDIR)/pbsa.MPI$(SFX)

libinstall: libinstall.$(INSTALLTYPE)

libinstall.serial: libpbsa.a
	/bin/mv libpbsa.a $(LIBDIR)

libinstall.parallel:
	@(if [ ! -f $(LIBDIR)/libpbsa.a ]; then \
		echo "Error: a serial version of libpbsa.a must be built before parallel run." ;\
		exit 2;\
	  fi ;\
	)
#---------------------------------------------------------------------------
configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\
		fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) and do NOT specify -mpi." ;\
		exit 2 ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is of type serial, not parallel !" ;\
			echo "  Rerun $(CONFIG_COMMAND) and specify an MPI implementation." ;\
			exit 2 ;\
		} ;\
	)

#---------------------------------------------------------------------------
pbsa$(SFX): $(OBJ) syslib netlib c9x-complex configured_serial
	$(FC) $(FFLAGS) $(AMBERFFLAGS) -o pbsa$(SFX) $(OBJ) \
	        ../lib/nxtsec.o ../lib/random.o \
	        $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

pbsa.MPI$(SFX): $(OBJ) syslib netlib c9x-complex configured_parallel
	$(FC) $(FFLAGS) $(AMBERFFLAGS) -o pbsa.MPI$(SFX) $(OBJ) \
	        ../lib/nxtsec.o ../lib/random.o \
	        $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

simplepbsa$(SFX): simplepbsa.o gopt.o libpbsa.a sfflib syslib netlib c9x-complex
	$(FC) $(FFLAGS) $(AMBERFFLAGS) -o simplepbsa$(SFX) simplepbsa.o gopt.o \
		libpbsa.a $(LIBDIR)/libsff.a ../lib/nxtsec.o $(FLIBSF) \
		$(LDFLAGS) $(AMBERLDFLAGS)
	/bin/mv simplepbsa$(SFX) $(BINDIR)

libpbsa.a: $(LIBPBSAOBJS) syslib
	-rm -f libpbsa.a
	$(AR) libpbsa.a $(LIBPBSAOBJS) ../lib/random.o
	$(RANLIB) libpbsa.a

libFpbsa.a: libFpbsa.$(INSTALLTYPE)

libFpbsa.serial: sanderdecomp.mod $(SANDERPBSAOBJS)
	$(AR) libFpbsa.a $(SANDERPBSAOBJS)
	$(RANLIB) libFpbsa.a
	/bin/mv libFpbsa.a $(LIBDIR)
	/bin/rm -f decomp.o decomp.mod

libFpbsa.parallel: $(SANDERPBSAOBJS)
	@(if [ ! -f $(LIBDIR)/libFpbsa.a ]; then \
		echo "Error: a serial version of libFpbsa.a must be built before parallel build." ;\
		exit 2;\
	  fi ;\
	)

sanderdecomp.mod:
	/bin/rm -f decomp.o decomp.mod
	/bin/cp ../../../src/sander/decomp.mod .

syslib:
	cd ../lib && $(MAKE) nxtsec.o random.o

sfflib:
	cd ../sff && $(MAKE) install

netlib: 
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../blas && $(MAKE) $(BLAS)
	cd ../arpack && $(MAKE) install

c9x-complex:
	@if test $(C9XCOMPLEX) != "skip"; then \
	   cd ../c9x-complex && $(MAKE) libmc.a; \
	 fi

pb_init.o: pb_init.f
	$(FPP) $(FPPFLAGS) $(AMBERFPPFLAGS) $< > _$<
	$(FC) -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ _$<

pb_init.SANDER.o: pb_init.f
	$(FPP) $(FPPFLAGS) $(AMBERFPPFLAGS) -DSANDER $< > _$<
	$(FC) -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ _$<

%.LIBPBSA.o: %.f
	$(FPP) $(FPPFLAGS) $(AMBERFPPFLAGS) -DLIBPBSA $< > _$<
	$(FC) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ _$<

%.SANDER.o: %.f
	$(FPP) $(FPPFLAGS) $(AMBERFPPFLAGS) -DSANDER $< > _$<
	$(FC) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ _$<

.f.o:   $<
	$(FPP) $(FPPFLAGS) $(AMBERFPPFLAGS) $< > _$<
	$(FC) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ _$<

clean:
	/bin/rm -f pbsa$(SFX) libpbsa.a
	/bin/rm -f *.o
	/bin/rm -f _*.f
	/bin/rm -f *.mod *__genmod.f90
	/bin/rm -f *.d

distclean: clean



#----------HEADER DEPENDENCIES (currently manually written)

myopen.o: \
                    copyright.h
getcoor.o: \
                    copyright.h\
                    files.h
egb.o: \
                    copyright.h\
                    md.h\
                    parms.h\
                    timer.h
ene.o: \
                    copyright.h\
                    md.h\
                    box.h\
                    parms.h\
                    flocntrl.h
pb_write.o: \
                    copyright.h\
                    box.h\
                    files.h\
                    md.h\
                    pb_def.h\
                    extra.h\
                    pb_force.o\
                    np_force.o\
                    decomp.o

pb_write.SANDER.o: \
                    copyright.h\
                    box.h\
                    files.h\
                    md.h\
                    pb_def.h\
                    extra.h\
                    pb_force.SANDER.o\
                    np_force.SANDER.o
force.o: \
                    copyright.h\
                    egb.o\
                    pb_constants.h\
                    timer.h\
                    md.h\
                    box.h\
                    memory.h\
                    parms.h\
                    files.h\
                    flocntrl.h\
                    pb_force.o\
                    np_force.o
locmem.o: \
                    copyright.h\
                    box.h\
                    memory.h\
                    md.h\
                    dynph.h
pb_read.o: \
                    copyright.h\
                    box.h\
                    pb_constants.h\
                    files.h\
                    md.h\
                    memory.h\
                    parms.h\
                    dynph.h\
                    pb_force.o\
                    np_force.o\
                    sa_driver.o\
                    decomp.o
pb_read.SANDER.o: \
                    copyright.h\
                    box.h\
                    pb_constants.h\
                    files.h\
                    md.h\
                    memory.h\
                    parms.h\
                    dynph.h\
                    pb_force.SANDER.o\
                    np_force.SANDER.o\
                    sa_driver.SANDER.o
pb_read.LIBPBSA.o: \
                    copyright.h\
                    box.h\
                    pb_constants.h\
                    timer.h\
                    files.h\
                    md.h\
                    memory.h\
                    parms.h\
                    dynph.h\
                    pb_force.o\
                    np_force.o
timer.o: \
                    copyright.h\
                    timer.h
sa_driver.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h\
                    memory.h\
                    timer.h\
                    extra.h\
                    timer.o
sa_driver.SANDER.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h\
                    memory.h\
                    timer.h\
                    extra.h\
                    timer.SANDER.o
np_force.o: \
                    copyright.h\
                    pb_def.h\
                    timer.h\
                    pb_md.h\
                    pb_constants.h\
                    pb_force.o\
		    sa_driver.o
np_force.SANDER.o: \
                    copyright.h\
                    pb_def.h\
                    pb_md.h\
                    pb_constants.h\
                    pb_force.SANDER.o\
		    sa_driver.SANDER.o
debug.o: \
                    copyright.h\
                    debug.h\
                    flocntrl.h\
                    md.h\
                    memory.h
decomp.o: \
                    copyright.h
rgroup.o: \
                    copyright.h
rfree.o: \
                    copyright.h
pb_force.o: \
                    copyright.h\
                    pb_def.h\
                    timer.h\
                    md.h\
                    pb_md.h\
                    box.h\
                    pb_constants.h\
                    pb_dbfrc_fld.h\
                    pb_dbfrc_fld2.h\
                    extra.h\
                    parallel.h\
		    sa_driver.o\
                    pb_fddrv.o\
                    pb_iimdrv.o\
                    timer.o\
                    decomp.o
pb_force.SANDER.o: \
                    copyright.h\
                    pb_def.h\
                    timer.h\
                    md.h\
                    pb_md.h\
                    box.h\
                    pb_constants.h\
                    pb_dbfrc_fld.h\
                    pb_dbfrc_fld2.h\
                    extra.h\
                    parallel.h\
		    sa_driver.SANDER.o\
                    pb_fddrv.o\
                    pb_iimdrv.o\
                    timer.SANDER.o\
                    decomp.o
pb_fddrv.o: \
                    copyright.h\
                    pb_def.h\
                    pb_constants.h\
                    pb_nlsolver.o\
                    pb_lsolver.o
pb_iimdrv.o: \
                    copyright.h\
                    pb_def.h\
                    pb_constants.h\
                    IIM.o
pb_exmol.o: \
                    copyright.h\
                    pb_def.h\
                    pb_constants.h\
                    md.h\
                    pb_force.o\
                    sa_driver.o\
                    timer.h\
                    timer.o
pb_exmol.SANDER.o: \
                    copyright.h\
                    pb_def.h\
                    pb_force.SANDER.o\
                    sa_driver.SANDER.o\
                    timer.h\
                    timer.SANDER.o
pb_mpfrc.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h
pb_mpfrc.SANDER.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h
pb_direct.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h\
                    pb_force.o\
                    decomp.o
pb_direct.SANDER.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h\
                    pb_force.SANDER.o
pb_list.o: \
                    copyright.h\
                    pb_def.h\
                    pb_constants.h\
                    extra.h\
                    pb_force.o\
                    sa_driver.o
pb_list.SANDER.o: \
                    copyright.h\
                    pb_def.h\
                    pb_constants.h\
                    extra.h\
                    pb_force.SANDER.o\
                    sa_driver.SANDER.o
pb_nlsolver.o: \
                    copyright.h\
                    pb_constants.h
pb_lsolver.o: \
                    copyright.h\
                    pb_constants.h
pb_init.o: \
                    copyright.h\
                    pb_def.h\
                    files.h\
                    md.h\
                    pb_md.h\
                    parms.h\
                    pb_constants.h\
                    sa_driver.o\
                    pb_force.o\
                    np_force.o
pb_init.SANDER.o: \
                    copyright.h\
                    pb_def.h\
                    files.h\
                    md.h\
                    pb_md.h\
                    parms.h\
                    pb_constants.h\
                    sa_driver.SANDER.o\
                    pb_force.SANDER.o\
                    np_force.SANDER.o
runmd.o: \
                    copyright.h\
                    md.h\
                    box.h\
                    memory.h\
                    timer.h\
                    pb_constants.h\
                    pb_md.h\
                    extra.h\
                    parms.h
runmin.o: \
                    copyright.h\
                    pb_constants.h\
                    md.h\
                    box.h\
                    files.h\
                    memory.h\
                    extra.h\
                    pb_md.h
pbsa.o: \
                    copyright.h \
                    files.h\
                    memory.h\
                    box.h\
                    md.h\
                    parms.h\
                    timer.h\
                    extra.h\
                    egb.o\
                    timer.o\
                    decomp.o
pbsaapi.o: \
                    pb_force.o\
		    np_force.o
pb_chkfrc.o: \
                    copyright.h\
                    pb_constants.h\
                    ../lib/random.h\
                    md.h\
                    pb_md.h\
                    memory.h\
                    pb_force.o
problem.o: \
                    pb_force.o
# DO NOT DELETE

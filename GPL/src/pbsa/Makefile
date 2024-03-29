#
#************************************************************************
#
include ../config.h

# PBSA Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

SRC=dsvdc.f  pbsa.f pb_read.f pb_write.f getcoor.f \
	runmd.f runmin.f force.f \
	pb_init.f sa_driver.f np_force.f pb_force.f pb_nlsolver.f \
	pb_exmol.f pb_mpfrc.f pb_direct.f pb_list.f timer.f \
	egb.f ene.f locmem.f myopen.f rdparm.f pb_lsolver.f \
	decomp.f rgroup.f rfree.f

OBJ=dsvdc.o  pbsa.o pb_read.o pb_write.o getcoor.o \
	runmd.o runmin.o force.o \
	pb_init.o sa_driver.o np_force.o pb_force.o pb_nlsolver.o \
	pb_exmol.o pb_mpfrc.o pb_direct.o pb_list.o timer.o \
	egb.o ene.o locmem.o myopen.o rdparm.o pb_lsolver.o \
	decomp.o rgroup.o rfree.o

LIBPBSAOBJS=interface.o pbsaapi.o \
        myopen.o pb_read.LIBPBSA.o pb_force.o sa_driver.o pb_nlsolver.o \
        pb_lsolver.o dsvdc.o timer.o pb_init.o pb_exmol.o np_force.o \
        pb_direct.o pb_list.o pb_write.o pb_mpfrc.o \
        decomp.o rgroup.o rfree.o

all: pbsa$(SFX)

install: pbsa$(SFX)
	/bin/mv pbsa$(SFX) $(BINDIR)

libinstall: libpbsa.a
	/bin/mv libpbsa.a $(LIBDIR)

uninstall:
	-rm -f $(BINDIR)/pbsa$(SFX) $(BINDIR)/simplepbsa$(SFX)

#---------------------------------------------------------------------------
pbsa$(SFX): $(OBJ) syslib netlib c9x-complex
	$(FC) -o pbsa$(SFX) $(OBJ) ../lib/nxtsec.o $(FLIBSF) $(LDFLAGS)

simplepbsa$(SFX): simplepbsa.o gopt.o libpbsa.a sfflib syslib netlib c9x-complex
	$(FC) -o simplepbsa$(SFX) simplepbsa.o gopt.o libpbsa.a \
	$(LIBDIR)/libsff.a ../lib/nxtsec.o $(FLIBSF) $(LDFLAGS)
	/bin/mv simplepbsa$(SFX) $(BINDIR)

libpbsa.a: $(LIBPBSAOBJS)
	-rm -f libpbsa.a
	$(AR) libpbsa.a $(LIBPBSAOBJS)
	$(RANLIB) libpbsa.a

syslib:
	cd ../lib && $(MAKE) nxtsec.o

sfflib:
	cd ../sff && $(MAKE) install

netlib: 
	cd ../blas && $(MAKE) $(BLAS)

c9x-complex:
	@if test $(C9XCOMPLEX) != "skip"; then \
	   cd ../c9x-complex && $(MAKE) libmc.a; \
	 fi

.f.o:   $<
	$(FPP) $< > _$<
	$(FC) -c $(FOPTFLAGS) -o $@ _$<

clean:
	/bin/rm -f pbsa$(SFX) 
	/bin/rm -f *.o
	/bin/rm -f _*.f
	/bin/rm -f *.mod *__genmod.f90
	/bin/rm -f *.d

distclean: clean

#----------HEADER DEPENDENCIES (constructed from ./makedepend)

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
                    extra.h\
                    pb_force.o\
                    np_force.o\
                    decomp.o
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
                    extra.h\
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
                    timer.h\
                    extra.h\
                    files.h\
                    md.h\
                    memory.h\
                    parms.h\
                    dynph.h\
                    pb_force.o\
                    np_force.o\
                    decomp.o\
                    rgroup.o
pb_read.LIBPBSA.o: \
                    copyright.h\
                    box.h\
                    pb_constants.h\
                    timer.h\
                    extra.h\
                    files.h\
                    md.h\
                    memory.h\
                    parms.h\
                    dynph.h\
                    pb_force.o\
                    np_force.o

pb_read.LIBPBSA.o: pb_read.f
	$(FPP) -DLIBPBSA pb_read.f > _pb_read.f
	$(FC) -c $(FOPTFLAGS) -o pb_read.LIBPBSA.o _pb_read.f

timer.o: \
                    copyright.h\
                    timer.h
sa_driver.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h
np_force.o: \
                    copyright.h\
                    pb_def.h\
                    timer.h\
                    pb_md.h\
                    pb_constants.h\
                    pb_force.o\
		    sa_driver.o
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
		    sa_driver.o\
                    pb_nlsolver.o\
                    pb_lsolver.o\
                    decomp.o
pb_exmol.o: \
                    copyright.h\
                    pb_def.h\
                    pb_constants.h
pb_mpfrc.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h
pb_direct.o: \
                    copyright.h\
                    pb_constants.h\
                    pb_def.h\
                    decomp.o
pb_list.o: \
                    copyright.h\
                    pb_def.h\
                    pb_constants.h\
                    pb_force.o
pb_nlsolver.o: \
                    copyright.h\
                    timer.h\
                    timer.o\
                    pb_constants.h
pb_lsolver.o: \
                    copyright.h\
                    timer.h\
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
runmd.o: \
                    copyright.h\
                    files.h\
                    md.h\
                    box.h\
                    memory.h\
                    extra.h\
                    timer.h\
                    pb_constants.h\
                    pb_md.h\
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
                    extra.h\
                    timer.h\
                    egb.o\
                    timer.o\
                    decomp.o
pbsaapi.o: \
                    pb_force.o\
		    np_force.o
# DO NOT DELETE

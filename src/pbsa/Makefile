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
	egb.f ene.f locmem.f myopen.f rdparm.f pb_lsolver.f

OBJ=dsvdc.o  pbsa.o pb_read.o pb_write.o getcoor.o \
	runmd.o runmin.o force.o \
	pb_init.o sa_driver.o np_force.o pb_force.o pb_nlsolver.o \
	pb_exmol.o pb_mpfrc.o pb_direct.o pb_list.o timer.o \
	egb.o ene.o locmem.o myopen.o rdparm.o pb_lsolver.o

all::	pbsa$(SFX)

install:	pbsa$(SFX)
	/bin/mv pbsa$(SFX) $(BINDIR)

uninstall:
	-rm -f $(BINDIR)/pbsa$(SFX)

#---------------------------------------------------------------------------
pbsa$(SFX):	$(OBJ) syslib netlib
	$(FC) -o pbsa$(SFX) $(OBJ) ../lib/sys2.a $(FLIBSF) $(LDFLAGS)

syslib:
	cd ../lib && $(MAKE) sys2.a

netlib: 
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../blas && $(MAKE) $(BLAS)

.f.o:   $<
	$(FPP) $< > _$<
	$(FC) -c $(FOPTFLAGS) -o $@ _$<

clean:
	/bin/rm -f pbsa$(SFX) 
	/bin/rm -f *.o
	/bin/rm -f _*.f
	/bin/rm -f *.mod *__genmod.f90
	/bin/rm -f *.d

#----------HEADER DEPENDENCIES (constructed from ./makedepend)

myopen.o: \
                    copyright.h
getcoor.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    files.h
egb.o: \
                    copyright.h\
                    ../include/assert.fh\
                    ../include/dprec.fh\
                    md.h\
                    parms.h\
                    timer.h\
                    constants.h
ene.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    md.h\
                    box.h\
                    parms.h\
                    flocntrl.h
assert.o: \
                    copyright.h\
                    ../include/assert.fh
pb_write.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    box.h\
                    files.h\
                    md.h\
                    extra.h\
                    pb_force.o\
                    np_force.o
findmask.o: \
                    ../include/dprec.fh
force.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    egb.o\
                    constants.h\
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
                    ../include/dprec.fh\
                    box.h\
                    memory.h\
                    md.h\
                    dynph.h
pb_read.o: \
                    copyright.h\
                    box.h\
                    constants.h\
                    timer.h\
                    extra.h\
                    files.h\
                    md.h\
                    memory.h\
                    parms.h\
                    dynph.h\
                    ../include/dprec.fh\
                    ../lib/random.h\
                    findmask.o\
                    pb_force.o\
                    np_force.o
timer.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    timer.h
sa_driver.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    constants.h\
                    pb_def.h
np_force.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    pb_def.h\
                    timer.h\
                    pb_md.h\
                    constants.h\
                    pb_force.o\
		    sa_driver.o
pb_force.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    pb_def.h\
                    timer.h\
                    md.h\
                    pb_md.h\
                    box.h\
                    constants.h\
                    pb_dbfrc_fld.h\
                    pb_dbfrc_fld2.h\
		    sa_driver.o\
                    pb_nlsolver.o\
                    pb_lsolver.o
pb_exmol.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    pb_def.h\
                    constants.h
pb_mpfrc.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    constants.h\
                    pb_def.h
pb_direct.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    constants.h\
                    pb_def.h
pb_list.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    pb_def.h\
                    constants.h\
                    pb_force.o
pb_nlsolver.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    timer.h\
                    timer.o\
                    constants.h
pb_lsolver.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    timer.h\
                    constants.h
pb_init.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    pb_def.h\
                    files.h\
                    md.h\
                    pb_md.h\
                    parms.h\
                    constants.h\
                    sa_driver.o\
                    pb_force.o\
                    np_force.o
runmd.o: \
                    copyright.h\
                    ../include/dprec.fh\
                    ../include/assert.fh\
                    files.h\
                    md.h\
                    box.h\
                    memory.h\
                    extra.h\
                    timer.h\
                    constants.h\
                    pb_md.h\
                    ../lib/random.h\
                    parms.h
runmin.o: \
                    copyright.h\
                    ../include/assert.fh\
                    ../include/dprec.fh\
                    constants.h\
                    md.h\
                    box.h\
                    files.h\
                    memory.h\
                    extra.h\
                    pb_md.h
pbsa.o: \
                    copyright.h \
                    ../include/dprec.fh\
                    ../include/assert.fh\
                    files.h\
                    memory.h\
                    box.h\
                    md.h\
                    parms.h\
                    extra.h\
                    timer.h\
                    egb.o\
                    timer.o
# DO NOT DELETE

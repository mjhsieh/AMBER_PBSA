#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
#                Regents of the University of California               **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
#
include ../config.h

# Sander Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

SRC=	\
	constants.f pbsa.f runmd.f runmin.f mdread.f getparm.f getcor.f \
	force.f is_init.f sa_driver.f np_force.f pb_force.f \
        pb_exmol.f pb_mpfrc.f pb_direct.f pb_list.f \
	egb.f ene.f decomp.f icosasurf.f findmask.f printe.f locmem.f \
	fastwt.f echoin.f parallel.f shake.f mdwrit.f minrit.f \
	set.f setmm.f dynlib.f mdfil.f nmlsrc.f \
	rstack.f istack.f rfree.f rgroup.f \
	random.f amopen.f debug.f degcnt.f \
	mexit.f new_time.f matinv.f assert.f veclib.f constantph.f \
        prm.c

OBJ=	\
	constants.o pbsa.o runmd.o runmin.o mdread.o getparm.o getcor.o \
	force.o is_init.o sa_driver.o np_force.o pb_force.o \
        pb_exmol.o pb_mpfrc.o pb_direct.o pb_list.o \
	egb.o ene.o decomp.o icosasurf.o findmask.o printe.o locmem.o \
	fastwt.o echoin.o parallel.o shake.o mdwrit.o minrit.o \
	set.o setmm.o dynlib.o mdfil.o nmlsrc.o \
	rstack.o istack.o rfree.o rgroup.o \
	random.o amopen.o debug.o degcnt.o \
	mexit.o new_time.o matinv.o assert.o veclib.o \
        prm.o


all::    pbsa$(SFX)

install:	pbsa$(SFX)
	/bin/mv pbsa$(SFX) $(BINDIR)

uninstall:
	-rm -f $(BINDIR)/pbsa$(SFX)

#---------------------------------------------------------------------------
pbsa$(SFX):	$(OBJ) syslib nxtsec netlib
	$(FC) -o pbsa$(SFX) $(OBJ) ../lib/sys.a ../lib/nxtsec.o $(FLIBSF) $(LDFLAGS) ../../lib/libsff.a

#-----------LIBS

syslib:
	cd ../lib && $(MAKE) sys.a

netlib:  
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../blas && $(MAKE) $(BLAS)

nxtsec:  
	cd ../lib && $(MAKE) nxtsec.o

.f.o:   $<
	$(FPP) $< > _$<
	$(FC) -c $(FFLAGS) -o $@ _$<

.c.o:	$<
	$(CC) -c $(CFLAGS) -o $@ $<

clean:
	/bin/rm -f pbsa$(SFX) 
	/bin/rm -f *.o
	/bin/rm -f _*.f
	/bin/rm -f *.mod
	/bin/rm -f *.d

#----------HEADER DEPENDENCIES (constructed from ./makedepend)

amopen.o: \
                    copyright.h
assert.o: \
                    copyright.h\
                    ../include/assert.fh
constantph.o: \
                    copyright.h\
                    dprec.h\
                    dynph.h\
                    random.h\
                    md.h\
                    files.h\
                    constants.h\
                    extra.h
debug.o: \
                    copyright.h\
                    dprec.h\
                    flocntrl.h\
                    debug.h\
                    md.h
decomp.o: \
                    copyright.h\
                    ../include/assert.fh\
                    dprec.h\
                    memory.h
degcnt.o: \
                    copyright.h\
                    dprec.h
dynlib.o: \
                    copyright.h\
                    dprec.h\
                    box.h\
                    files.h\
                    md.h\
                    extra.h\
                    parallel.h
echoin.o: \
                    copyright.h
egb.o: \
                    copyright.h\
                    ../include/assert.fh\
                    dprec.h\
                    icosasurf.o\
                    decomp.o\
                    parallel.h\
                    md.h\
                    parms.h\
                    def_time.h\
                    constants.h\
                    gbsa.h
ene.o: \
                    copyright.h\
                    dprec.h\
                    decomp.o\
                    parallel.h\
                    md.h\
                    box.h\
                    parms.h\
                    flocntrl.h
fastwt.o: \
                    copyright.h\
                    dprec.h\
                    parallel.h\
                    md.h\
                    constants.h
findmask.o: \
                    dprec.h
force.o: \
                    copyright.h\
                    dprec.h\
                    egb.o\
                    pb_force.o\
                    np_force.o\
                    constants.h\
                    def_time.h\
                    parallel.h\
                    md.h\
                    box.h\
                    memory.h\
                    parms.h\
                    files.h\
                    extra.h\
                    flocntrl.h
getcor.o: \
                    dprec.h\
                    copyright.h\
                    files.h
icosasurf.o: \
                    copyright.h\
                    ../include/assert.fh\
                    dprec.h\
                    decomp.o
istack.o: \
                    copyright.h\
                    memory.h\
                    istack.h
locmem.o: \
                    copyright.h\
                    dprec.h\
                    box.h\
                    memory.h\
                    md.h\
                    dynph.h\
                    parallel.h
matinv.o: \
                    copyright.h\
                    dprec.h
mdfil.o: \
                    copyright.h\
                    files.h\
                    parallel.h
mdread.o: \
                    copyright.h\
                    dprec.h\
                    box.h\
                    constants.h\
                    def_time.h\
                    extra.h\
                    files.h\
                    md.h\
                    memory.h\
                    parms.h\
                    dynph.h\
                    random.h\
                    decomp.o\
                    findmask.o
mdwrit.o: \
                    copyright.h\
                    dprec.h\
                    files.h
mexit.o: \
                    copyright.h\
                    parallel.h
minrit.o: \
                    copyright.h\
                    dprec.h\
                    md.h\
                    files.h\
                    box.h
new_time.o: \
                    copyright.h\
                    dprec.h\
                    def_time.h\
                    new_time.h\
                    extra.h\
                    istack.h\
                    rstack.h\
                    parallel.h
nmlsrc.o: \
                    copyright.h
parallel.o: \
                    copyright.h\
                    dprec.h\
                    parallel.h\
                    extra.h\
                    md.h\
                    memory.h\
                    parms.h\
                    box.h\
                    files.h\
                    nonper.h\
                    flocntrl.h\
                    debug.h\
                    rstack.h\
                    istack.h\
                    new_time.h
sa_driver.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    constants.o\
                    is_def.h
np_force.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    def_time.h\
                    pb_md.h\
                    constants.o\
                    pb_force.o\
		    sa_driver.o
pb_force.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    def_time.h\
                    md.h\
                    pb_md.h\
                    box.h\
                    constants.o\
		    sa_driver.o
pb_exmol.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    constants.o\
                    pb_force.o
pb_mpfrc.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    constants.o\
                    is_def.h
pb_direct.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    constants.o\
                    is_def.h
pb_list.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    constants.o\
                    pb_force.o
is_init.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    files.h\
                    md.h\
                    pb_md.h\
                    parms.h\
                    constants.o\
                    sa_driver.o\
                    np_force.o\
                    pb_force.o
printe.o: \
                    copyright.h\
                    dprec.h\
                    decomp.o\
                    box.h\
                    extra.h\
                    files.h\
                    md.h\
                    memory.h\
                    constants.h
random.o: \
                    copyright.h\
                    dprec.h\
                    random.h
getparm.o: \
                    copyright.h\
                    dprec.h\
                    md.h\
                    memory.h\
                    parms.h\
                    files.h\
                    box.h
rfree.o: \
                    copyright.h\
                    dprec.h
rgroup.o: \
                    copyright.h\
                    dprec.h
rstack.o: \
                    copyright.h\
                    memory.h\
                    rstack.h
runmd.o: \
                    copyright.h\
                    dprec.h\
                    ../include/assert.fh\
                    decomp.o\
                    fastwt.o\
                    parallel.h\
                    files.h\
                    md.h\
                    box.h\
                    memory.h\
                    extra.h\
                    def_time.h\
                    constants.h\
                    pb_md.h\
                    random.h\
                    parms.h
runmin.o: \
                    copyright.h\
                    ../include/assert.fh\
                    dprec.h\
                    fastwt.o\
                    parallel.h\
                    constants.h\
                    md.h\
                    box.h\
                    files.h\
                    memory.h\
                    extra.h\
                    pb_md.h
pbsa.o: \
                    copyright.h\
                    dprec.h\
                    ../include/assert.fh\
                    files.h\
                    parallel.h\
                    egb.o\
                    decomp.o\
                    fastwt.o\
                    memory.h\
                    box.h\
                    md.h\
                    parms.h\
                    extra.h\
                    nonper.h\
                    def_time.h
set.o: \
                    copyright.h\
                    dprec.h\
                    files.h\
                    md.h\
                    parms.h\
                    memory.h\
                    parallel.h\
                    extra.h
setmm.o: \
                    copyright.h\
                    dprec.h\
                    memory.h\
                    box.h
shake.o: \
                    copyright.h\
                    dprec.h\
                    parallel.h\
                    md.h\
                    extra.h
veclib.o: \
                    copyright.h\
                    dprec.h\
                    constants.h
prm.o: \
                    is_copyrightc.h\
                    prm.h

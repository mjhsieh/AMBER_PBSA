include ../config.h
SHELL=/bin/sh
SRCDIR=..
VPATH=

.f.o:
	$(FC) -c $(FOPTFLAGS) -o $@ $<

OBJS= lsame.o

install: $(OBJS)
	$(AR) lapack.a $(OBJS)
	$(RANLIB) lapack.a
	mv lapack.a $(LIBDIR)

skip:
	@echo "skipping compilation of LAPACK"

uninstall:
	rm -f $(LIBDIR)/lapack.a

clean:
	-/bin/rm -f *.o 

dlamch.o:  dlamch.f
	$(FC) -c $(FFLAGS) -o $@ $<

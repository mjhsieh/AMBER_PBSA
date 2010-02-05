include ../config.h
SHELL=/bin/sh
SRCDIR=..
VPATH=

.f.o:
	$(FPP) $< > _$<
	$(FC) -c $(FOPTFLAGS) -o $@ _$<

OBJS= drotg.o daxpy.o dswap.o drot.o dnrm2.o ddot.o dscal.o dgemm.o xerbla.o \
      dswap.o

install: $(OBJS)
	$(AR) blas.a $(OBJS)
	$(RANLIB) blas.a
	mv blas.a $(LIBDIR)

skip:
	@echo "skipping compilation of BLAS"

clean:
	-rm -f *.o _*.f

uninstall:
	-rm -f $(LIBDIR)/blas.a


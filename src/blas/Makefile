include ../config.h
SHELL=/bin/sh
SRCDIR=..
VPATH=

.f.o:
	$(FPP) $< > _$<
	$(FC) -c $(FOPTFLAGS) -o $@ _$<

OBJS=    dasum.o daxpy.o dcabs1.o dcopy.o ddot.o \
         dgbmv.o dgemm.o dgemv.o dger.o dnrm2.o drot.o \
         drotg.o dsbmv.o dscal.o dspmv.o dspr.o dspr2.o \
         dswap.o dsymm.o dsymv.o dsyr.o dsyr2.o dsyr2k.o \
         dsyrk.o dtbmv.o dtbsv.o dtpmv.o dtpsv.o dtrmm.o \
         dtrmv.o dtrsm.o dtrsv.o dummy.o dzasum.o dznrm2.o \
         idamax.o isamax.o lsame.o xerbla.o \
         scopy.o sswap.o sscal.o sasum.o snrm2.o sdot.o ssyr2k.o \
         saxpy.o sgemv.o ssyr2.o ssymv.o sgemm.o strmm.o strmv.o sger.o \
         zher2k.o zswap.o zgemv.o zhemv.o zscal.o zdotc.o zaxpy.o \
         zher2.o zdscal.o ztrmv.o zcopy.o ztrmm.o zgemm.o zgerc.o

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


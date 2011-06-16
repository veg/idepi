MAKE = make
CC = gcc-4.6
CFLAGS = -O3 -fno-common -ffast-math -fPIC -DUSE_SSE2
CONTRIB = contrib
PYDIR   = idepi
SUBDIRS = $(wildcard $(CONTRIB)/hmmer*)

compf: libcompf.dylib

libcompf.dylib: idepi/_compf.c
	$(CC) $(CFLAGS) -Icontrib/sse_mathfun -c -o _compf.o $<
	$(CC) -shared -dylib -o $@ _compf.o
	@-rm _compf.o

all:
	@$(foreach var, $(SUBDIRS), make -C $(var) all;)

clean:
	@-$(foreach var, $(SUBDIRS), make -C $(var) clean;)
	@-rm $(PYDIR)/*.pyc
	@-rm libcompf.dylib

distclean: clean
	@-$(foreach var, $(SUBDIRS), make -C $(var) distclean;)
	@-rm $(PYDIR)/*.pyc

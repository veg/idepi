MAKE = make

CONTRIB = contrib
PYDIR   = idepi
SUBDIRS = $(wildcard $(CONTRIB)/hmmer*)

all:
	@$(foreach var, $(SUBDIRS), make -C $(var) all;)

clean:
	@-$(foreach var, $(SUBDIRS), make -C $(var) clean;)
	@-rm $(PYDIR)/*.pyc

distclean: clean
	@-$(foreach var, $(SUBDIRS), make -C $(var) distclean;)
	@-rm $(PYDIR)/*.pyc

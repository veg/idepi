MAKE = make
CONTRIB = contrib
PYDIR   = idepi
SUBDIRS = $(wildcard $(CONTRIB)/hmmer*)

all:
	@$(foreach var, $(SUBDIRS), $(MAKE) -C $(var) all;)

clean:
	@-$(foreach var, $(SUBDIRS), $(MAKE) -C $(var) clean;)
	@-rm $(PYDIR)/*.pyc

distclean: clean
	@-$(foreach var, $(SUBDIRS), $(MAKE) -C $(var) distclean;)
	@-rm $(PYDIR)/*.pyc

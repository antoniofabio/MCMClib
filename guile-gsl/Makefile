PREFIX:=/usr
LIBDIR:=$(PREFIX)/lib

MODULES := gsl
MODULES_LIBS := $(MODULES:%=libguile%.so)
SWIG_TARGETS := $(MODULES:%=%_wrap.c) $(MODULES:%=swig/%.scm)

GUILE_LIBRARY_PATH:= $(shell guile -c "(display (%package-data-dir))")

SWIG_FLAGS := -guile -scm -nodefaultctor -package swig -Linkage passive -scmstub
GUILE_LDFLAGS := $(LDFLAGS) -lguile -lgsl -lgslcblas -lm
GUILE_CFLAGS := $(CFLAGS) -I./src

all: swig-targets $(MODULES_LIBS)

swig-targets: $(SWIG_TARGETS)

gsl_wrap.c: vector.i matrix.i blas.i rng.i qrng.i sf.i mode.i permutation.i \
	linalg.i

%_wrap.c swig/%.scm: %.i
	swig $(SWIG_FLAGS) -o $@ $<

libguile%.so: %_wrap.c
	$(CC) -shared $(GUILE_CFLAGS) -fPIC $^ $(GUILE_LDFLAGS) -o $@

install: all
	mkdir -p $(GUILE_LIBRARY_PATH)/swig
	install -p -m 0644 -t $(GUILE_LIBRARY_PATH)/swig swig/*
	mkdir -p $(LIBDIR)
	install -p -m 0755 -t $(LIBDIR) $(MODULES_LIBS)
	mkdir -p $(GUILE_LIBRARY_PATH)/srfi
	install -p -m 0644 -t $(GUILE_LIBRARY_PATH)/srfi srfi/*

clean:
	@rm -f *.o $(MODULES_LIBS)

clean-all: clean
	@rm -f $(SWIG_TARGETS)

VERSION := $(shell git describe --abbrev=5)
TARBALL_NAME := guile-gsl_$(VERSION)

distrib: swig-targets
	git archive --format=tar --prefix=$(TARBALL_NAME)/ HEAD > $(TARBALL_NAME).tar
	tar xf $(TARBALL_NAME).tar
	unlink $(TARBALL_NAME).tar
	cp gsl_wrap.c $(TARBALL_NAME)
	cp swig/gsl.scm $(TARBALL_NAME)/swig
	rm -f $(TARBALL_NAME).tar.gz
	tar czf $(TARBALL_NAME).tar.gz $(TARBALL_NAME)
	rm -rf $(TARBALL_NAME)

.PHONY : all clean clean-all install config swig-targets distrib

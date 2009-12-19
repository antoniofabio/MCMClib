PREFIX:=/usr
LIBDIR:=$(PREFIX)/lib

MODULES := gsl
MODULES_LIBS := $(MODULES:%=libguile%.so)
MODULES_TARGETS := $(MODULES_LIBS) $(MODULES:%=%_wrap.c) $(MODULES:%=swig/%.scm)

GUILE_LIBRARY_PATH:= $(shell guile -c "(display (%package-data-dir))")

SWIG_FLAGS := -guile -scm -nodefaultctor -package swig -Linkage passive -scmstub
GUILE_LDFLAGS := $(LDFLAGS) -lguile -lgsl -lgslcblas -lm
GUILE_CFLAGS := $(CFLAGS) -I./src

all: $(MODULES_TARGETS)

gsl_wrap.c: vector.i matrix.i blas.i rng.i qrng.i sf.i mode.i permutation.i \
	linalg.i

%_wrap.c: %.i
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
	@rm -f $(MODULES_TARGETS)

.PHONY : all clean clean-all install config

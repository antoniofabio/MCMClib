CFLAGS+=  --std=gnu99 -O0 -g -Wall -fpic
LDFLAGS+= -lgsl -lgslcblas -lm
PREFIX:=/usr
LIBDIR:=$(PREFIX)/lib
INCLUDEDIR:=$(PREFIX)/include

MODULES := src test examples guile

SRC :=
TOCLEAN :=
include $(MODULES:%=%/module.mk)

.PHONY : all lib test examples guile clean doc

all: lib examples

lib: src/libmcmclib.a src/libmcmclib.so

examples: $(EXAMPLES_BIN)

doc:
	cd doc; sh makeDoxyfile.sh; doxygen

clean:
	@rm -rf *~ $(TOCLEAN) doc/html doc/latex

VERSION:=$(shell git describe --abbrev=5)
TARBALL_NAME:=MCMClib_$(VERSION)

distrib: $(SRC) doc
	tar cf $(TARBALL_NAME).tar $(SRC) doc Makefile $(MODULES:%=%/module.mk)\
		ChangeLog README COPYING
	mkdir $(TARBALL_NAME)
	cd $(TARBALL_NAME) && tar xf ../$(TARBALL_NAME).tar
	rm $(TARBALL_NAME).tar
	tar czf $(TARBALL_NAME).tar.gz $(TARBALL_NAME)
	rm -rf $(TARBALL_NAME)

install: all
	mkdir -p $(LIBDIR)
	install -p -m 0755 src/libmcmclib.so $(LIBDIR)
	install -p -m 0644 src/libmcmclib.a $(LIBDIR)
	mkdir -p $(INCLUDEDIR)
	install -p -m 0644 -t $(INCLUDEDIR) src/*.h

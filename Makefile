CFLAGS+=  -ansi -pedantic -std=c99 -Wall -W  -Wconversion -Wshadow -Wpointer-arith \
	-Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums \
	-fno-common -Dinline= -g -O2 -fPIC
LDFLAGS+= -lgsl -lgslcblas -lm
PREFIX:=/usr
LIBDIR:=$(PREFIX)/lib
INCLUDEDIR:=$(PREFIX)/include

MODULES := src test guile test-guile R

SRC :=
TOCLEAN :=
include $(MODULES:%=%/module.mk)

.PHONY : all lib test guile clean doc test-guile

all: lib

lib: src/libmcmclib.a src/libmcmclib.so

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

check-syntax:
	$(CC) $(CFLAGS) -Isrc --std=c99 -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)

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
	@rm -rf *~ $(SRC:%.c=%.o) $(SRC:%=%~) $(TOCLEAN) doc/html doc/latex

distrib:
	./distrib.sh

install: all
	mkdir -p $(LIBDIR)
	install -p -m 0755 src/libmcmclib.so $(LIBDIR)
	install -p -m 0644 src/libmcmclib.a $(LIBDIR)
	mkdir -p $(INCLUDEDIR)
	install -p -m 0644 -t $(INCLUDEDIR) src/*.h

CFLAGS+=  --std=gnu99 -O0 -g -Wall -fpic
LDFLAGS+= -lgsl -lgslcblas -lm

MODULES := src test examples

SRC :=
TOCLEAN :=
include $(MODULES:%=%/module.mk)

.PHONY : all lib test examples clean doc

all: lib examples

lib: src/libmcmclib.a src/libmcmclib.so

examples: $(EXAMPLES_BIN)

doc:
	cd doc; doxygen

clean:
	@rm -rf *~ $(SRC:%.c=%.o) $(SRC:%=%~) $(TOCLEAN) doc/html doc/latex

distrib:
	./distrib.sh

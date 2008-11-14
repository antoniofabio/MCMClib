CFLAGS+=  --std=gnu99 -O0 -g -Wall -fpic
LDFLAGS+= -lgsl -lgslcblas -lm

MODULES := src examples

SRC :=
TOCLEAN :=
include $(MODULES:%=%/module.mk)

.PHONY : all lib examples clean

all: lib examples

lib: src/libmcmclib.a src/libmcmclib.so

examples: $(EXAMPLES_BIN)

clean:
	rm -rf $(SRC:%.c=%.o) $(SRC:%=%~) $(TOCLEAN)

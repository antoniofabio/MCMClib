CFLAGS+=  --std=gnu99 -O0 -g -Wall
LDFLAGS+= -lgsl -lgslcblas -lm

SOURCES = $(wildcard *.c)
OBJS = $(addsuffix .o, $(basename $(SOURCES)))

all: libmcmclib.a libmcmclib.so

libmcmclib.a: $(OBJS)
	ar rcs $@ $<

libmcmclib.so: $(OBJS)
	ld -shared $< $(LDFLAGS) -o $@

clean:
	rm -rf *.o *~ libmcmclib.a libmcmclib.so

CFLAGS+=  --std=gnu99 -O0 -g -Wall -fpic
LDFLAGS+= -lgsl -lgslcblas -lm

SOURCES = $(wildcard *.c)
OBJS = $(addsuffix .o, $(basename $(SOURCES)))

all: libmcmclib.a libmcmclib.so

libmcmclib.a: $(OBJS)
	ar rcs $@ $(OBJS)

libmcmclib.so: $(OBJS)
	$(CC) -shared $(OBJS) $(LDFLAGS) -o $@

clean:
	rm -rf *.o *~ libmcmclib.a libmcmclib.so

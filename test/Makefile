SRCDIR= ../src
CFLAGS+=  --std=gnu99 -O0 -g -Wall -I$(SRCDIR)
LDFLAGS+= $(SRCDIR)/libmcmclib.a -lgsl -lgslcblas -lm
TEST_LIST= t1 t2 tmixem

.PHONY: all test clean

all: $(TEST_LIST)

test: $(TEST_LIST)
	for i in $(TEST_LIST); do ./$$i; done

%: %.c $(SRCDIR)/libmcmclib.a
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

clean:
	rm -rf *.o *~ $(TEST_LIST)

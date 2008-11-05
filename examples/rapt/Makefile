SRCDIR= ../../src
CFLAGS+=  --std=gnu99 -O0 -g -Wall -I$(SRCDIR)
LDFLAGS+= $(SRCDIR)/libmcmclib.a -lgsl -lgslcblas -lm
EX_LIST= ex1 ex2 ex3

all: $(EX_LIST)

run: $(EX_LIST)
	./ex1
	./ex2

%: %.c $(SRCDIR)/libmcmclib.a
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

clean:
	rm -rf *.o *~ $(EX_LIST)

clean_all: clean
	rm -rf *.csv

SRCDIR= ../../src
CFLAGS=  --std=gnu99 -O0 -g -Wall -I$(SRCDIR)
LDFLAGS= -lgsl -lgslcblas -lm
EX_LIST= ex1

ex1: ex1.o $(SRCDIR)/libmcmclib.a

clean:
	rm -rf *.o $(EX_LIST)

clean_all: clean
	rm -rf *.csv

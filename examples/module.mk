EXES := spatial/example1 wrapped/example1 wrapped/regression\
	mcar/example1 poisson/example1 mcar/example2 mcar/example3\
	mcar/example4

EXAMPLES_BIN := $(EXES:%=examples/%)
TOCLEAN += $(EXAMPLES_BIN)
EXAMPLES_LDFLAGS := src/libmcmclib.a $(LDFLAGS)
EXAMPLES_CFLAGS := $(CFLAGS) -I./src

%.o: %.c
	$(CC) -c $(EXAMPLES_CFLAGS) $^ -o $@

$(EXAMPLES_BIN): %: %.c src/libmcmclib.a
	$(CC) $(EXAMPLES_CFLAGS) $^ $(EXAMPLES_LDFLAGS) -o $@


EXES := ex1 ex2 ex3 ex4

EXAMPLES_BIN := $(EXES:%=examples/rapt/%)
TOCLEAN += $(EXAMPLES_BIN)
EXAMPLES_LDFLAGS := src/libmcmclib.a $(LDFLAGS)
EXAMPLES_CFLAGS := $(CFLAGS) -I./src

$(EXAMPLES_BIN): %: %.c src/libmcmclib.a
	$(CC) $(EXAMPLES_CFLAGS) $< $(EXAMPLES_LDFLAGS) -o $@

examples/rapt/ex3: examples/rapt/ex3_target_distrib.c
examples/rapt/ex4: examples/rapt/ex4_target_distrib.c

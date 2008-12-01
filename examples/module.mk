EXES := example1 example2 example3
RAPT_EXES := ex2 

EXAMPLES_RAPT_BIN := $(RAPT_EXES:%=examples/rapt/%)
EXAMPLES_BIN := $(EXES:%=examples/%) $(EXAMPLES_RAPT_BIN)
TOCLEAN += $(EXAMPLES_BIN)
EXAMPLES_LDFLAGS := src/libmcmclib.a $(LDFLAGS)
EXAMPLES_CFLAGS := $(CFLAGS) -I./src

$(EXAMPLES_BIN): %: %.c src/libmcmclib.a
	$(CC) $(EXAMPLES_CFLAGS) $< $(EXAMPLES_LDFLAGS) -o $@

examples/rapt/ex3: examples/rapt/ex3_target_distrib.c
examples/rapt/ex4: examples/rapt/ex4_target_distrib.c

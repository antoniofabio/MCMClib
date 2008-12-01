EXES := example1 example2 example3
RAPT_EXES := ex2 olem

EXAMPLES_RAPT_BIN := $(RAPT_EXES:%=examples/rapt/%)
EXAMPLES_BIN := $(EXES:%=examples/%) $(EXAMPLES_RAPT_BIN)
TOCLEAN += $(EXAMPLES_BIN)
EXAMPLES_LDFLAGS := src/libmcmclib.a $(LDFLAGS)
EXAMPLES_CFLAGS := $(CFLAGS) -I./src

%.o: %.c
	$(CC) -c $(EXAMPLES_CFLAGS) $^ -o $@

$(EXAMPLES_BIN): %: %.c src/libmcmclib.a
	$(CC) $(EXAMPLES_CFLAGS) $^ $(EXAMPLES_LDFLAGS) -o $@

examples/rapt/olem: examples/rapt/mixnorm_target.o

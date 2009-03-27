EXES := spatial/example1

EXAMPLES_BIN := $(EXES:%=examples/%)
TOCLEAN += $(EXAMPLES_BIN)
EXAMPLES_LDFLAGS := src/libmcmclib.a $(LDFLAGS)
EXAMPLES_CFLAGS := $(CFLAGS) -I./src

%.o: %.c
	$(CC) -c $(EXAMPLES_CFLAGS) $^ -o $@

$(EXAMPLES_BIN): %: %.c src/libmcmclib.a
	$(CC) $(EXAMPLES_CFLAGS) $^ $(EXAMPLES_LDFLAGS) -o $@


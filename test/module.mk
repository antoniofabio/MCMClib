TEST_CFLAGS := $(CFLAGS) -I./src
TEST_LDFLAGS:= src/libmcmclib.a $(LDFLAGS)
TEST_names:= t1 t2 trapt tmixem tmixem_rec tmixem_online tolemrapt
TEST_BIN:= $(TEST_names:%=test/%)

TOCLEAN += $(TEST_BIN)

test: $(TEST_BIN)
	cd test; for i in $(TEST_names); do ./$$i; done

$(TEST_BIN): %: %.c src/libmcmclib.a
	$(CC) $(TEST_CFLAGS) $< $(TEST_LDFLAGS) -o $@

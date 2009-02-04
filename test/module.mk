TEST_CFLAGS := $(CFLAGS) -I./src
TEST_LDFLAGS:= src/libmcmclib.a $(LDFLAGS)
TEST_names:= t1 t2 trecursive_variance trapt tmixem tmixem_rec tmixem_online tolemrapt\
	trapt_inca tam_inca tolemrapt_inca
TEST_targets:=$(TEST_names:%=test_%)
TEST_BIN:= $(TEST_names:%=test/%)

TOCLEAN += $(TEST_BIN)

test: $(TEST_targets)

$(TEST_targets): test_%: test/%
	@cd test; echo -n "$<: "; if ../$<; then echo OK; else echo FAIL; fi

$(TEST_BIN): %: %.c src/libmcmclib.a
	$(CC) $(TEST_CFLAGS) $< $(TEST_LDFLAGS) -o $@

TEST_CFLAGS := $(CFLAGS) -I./src
TEST_LDFLAGS:= src/libmcmclib.a $(LDFLAGS)
TEST_names:= t1 t2 trecursive_variance tmixem tmixem_rec tmixem_online \
	tmh_q tmh tmh2 tgauss_rw tgauss_mrw \
	tamh tamh2 trapt_q trapt traptor tgauss_am \
	tinca trapt_inca tam_inca tinca_raptor \
	tspatial tscam tmcar_tilde twishart tmcar_model tiwishart \
	tgivens tmonitor

TEST_targets:=$(TEST_names:%=test_%)
TEST_BIN:= $(TEST_names:%=test/%)

TOCLEAN += $(TEST_BIN)

test: $(TEST_targets)

$(TEST_targets): test_%: test/%
	@cd test; echo -n "$<: "; if ../$<; then echo OK; else echo FAIL; fi

$(TEST_BIN): %: %.c src/libmcmclib.a
	$(CC) $(TEST_CFLAGS) $< $(TEST_LDFLAGS) -o $@

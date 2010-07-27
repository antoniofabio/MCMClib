# TEST_names:= tinca tmcar_tilde tgivens thier tpois_model
TEST_names := t1 t2 \
	twishart tiwishart tmcar_model \
	trecursive_variance tmixem tmixem_rec tmixem_online tmixem_online2 \
	tmh_q tmh tmh2 tmh3 tgauss_rw tgauss_mrw \
	tamh tamh2 trapt_q trapt traptor tgauss_am tat7 \
	tmonitor tmonitor2
TEST_CFLAGS := $(CFLAGS) -I./src
TEST_OBJ:= $(TEST_names:%=test/%.o)
TEST_LDFLAGS:= test/CuTest.o $(TEST_OBJ) src/libmcmclib.a $(LDFLAGS)
SRC += $(wildcard test/*.c)
SRC += $(wildcard test/*.h)
SRC += $(wildcard test/*.dat)
SRC += $(wildcard test/*.check)

TEST_BIN:= test/AllTests

TOCLEAN += $(TEST_BIN)

test: test/AllTests
	cd test; ./AllTests

$(TEST_OBJ):%.o: %.c
	$(CC) -c $< $(TEST_CFLAGS) -o $@

test/AllTests.c: $(TEST_names:%=test/%.c)
	cd test; ./make-tests.sh > AllTests.c

test/AllTests: test/AllTests.c src/libmcmclib.a test/CuTest.o $(TEST_OBJ)
	$(CC) $(TEST_CFLAGS) $< $(TEST_LDFLAGS) -o $@

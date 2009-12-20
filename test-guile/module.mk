SRC += $(wildcard test-guile/*.test)
SRC += $(wildcard test-guile/*.dat)
SRC += $(wildcard test-guile/*.scm)

test-guile:
	GUILE_LOAD_PATH=. ./test-guile/guile-test --test-suite test-guile\
		mh_q.test mh.test amh.test

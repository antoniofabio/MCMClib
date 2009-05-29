GUILE_MODULES := gsl mcmclib
GUILE_MODULES_TARGETS := $(GUILE_MODULES:%=guile/swig/libguile%.so)\
	$(GUILE_MODULES:%=guile/%_wrap.c)

SWIG_FLAGS := -guile -scm -nodefaultctor -package swig -Linkage passive -scmstub
GUILE_LDFLAGS := $(LDFLAGS) -lguile
GUILE_CFLAGS := $(CFLAGS) -I./src

TOCLEAN += $(GUILE_MODULES_TARGETS) $(GUILE_MODULES:%=guile/%_wrap.c)

guile: $(GUILE_MODULES_TARGETS)

guile/%_wrap.c: guile/%.i
	swig $(SWIG_FLAGS) -o $@ $<

guile/swig/libguile%.so: guile/%_wrap.c src/libmcmclib.a
	$(CC) -shared $(GUILE_CFLAGS) $(GUILE_LDFLAGS) -fPIC $< -o $@

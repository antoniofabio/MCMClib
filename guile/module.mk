GUILE_MODULES := mcmclib
GUILE_MODULES_LIBS := $(GUILE_MODULES:%=guile/swig/libguile%.so)
GUILE_MODULES_TARGETS := $(GUILE_MOULES_LIBS) $(GUILE_MODULES:%=guile/%_wrap.c)
GUILE_LIBRARY_PATH:= $(shell guile -c "(display (%package-data-dir))")

SWIG_FLAGS := -guile -scm -nodefaultctor -package swig -Linkage passive -scmstub
GUILE_LDFLAGS := $(LDFLAGS) -lguile
GUILE_CFLAGS := $(CFLAGS) -I./src

TOCLEAN += $(GUILE_MODULES:%=guile/swig/libguile%.so)

guile: $(GUILE_MODULES_TARGETS)

guile/mcmclib_wrap.c: guile/monitor.i guile/amh.i guile/distrfuns.i guile/mixem.i guile/at7.i
guile/%_wrap.c: guile/%.i
	swig $(SWIG_FLAGS) -o $@ $<

guile/swig/libguile%.so: guile/%_wrap.c src/libmcmclib.a
	$(CC) -shared $(GUILE_CFLAGS) -fPIC $^ $(GUILE_LDFLAGS) -o $@

install-guile: guile
	mkdir -p $(GUILE_LIBRARY_PATH)/swig
	install -p -m 0644 -t $(GUILE_LIBRARY_PATH)/swig guile/swig/*
	mkdir -p $(LIBDIR)
	install -p -m 0755 -t $(LIBDIR) $(GUILE_MODULES_LIBS)

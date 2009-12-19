LIB_SRC := $(wildcard src/*.c)
LIB_OBJ := $(LIB_SRC:%.c=%.o)

SRC += $(LIB_SRC) $(wildcard src/*.h)
TOCLEAN += src/libmcmclib.a src/libmcmclib.so

src/libmcmclib.a: $(LIB_OBJ)
	ar rcs $@ $(LIB_OBJ)

src/libmcmclib.so: $(LIB_OBJ)
	$(CC) -shared $(LIB_OBJ) $(LDFLAGS) -o $@

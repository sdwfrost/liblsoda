include ../Make.inc

SRCS =  \
cfode.c \
cfode_static.c \
common.c \
corfailure.c \
correction.c \
daxpy.c \
ddot.c \
dgefa.c \
dgesl.c \
dscal.c \
fnorm.c \
idamax.c \
intdy.c \
lsoda.c \
methodswitch.c \
orderswitch.c \
prja.c \
scaleh.c \
solsy.c \
stoda.c \
vmnorm.c \
strdup_printf.c

OBJS = $(SRCS:%.c=%.o)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS-add) $(CFLAGS) $(fPIC) -DNDEBUG -c $< -o $@

release debug: liblsoda.$(SHLIB_EXT)

liblsoda.$(SHLIB_EXT): $(OBJS)
	rm -rf $@
	$(CC) $(LDFLAGS) $(LDFLAGS-add) -shared -o $@ $^

clean:
	rm -f *.o *.$(SHLIB_EXT)

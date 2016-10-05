SRC =  \
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
main.c \
methodswitch.c \
orderswitch.c \
prja.c \
scaleh.c \
solsy.c \
stoda.c \
vmnorm.c \
strdup_printf.c

all: main.o liblsoda.a
	gcc -std=c99 -g -pg -O0 -o a.out $^ -lm -lgfortran
printcf: printcf.o liblsoda.a
	gcc -std=c99 -g -O0 -o printcf $^ -lm -lgfortran
test2: test2.o liblsoda.a
	gcc -std=c99 -g -pg -O0 -o test2 $^ -lm -lgfortran

clean:
	rm *.[oa]
liblsoda.a: $(SRC:.c=.o)
	ar cur $@ $^

.c.o:
	gcc -DCFODE_STATIC -std=c99 -g -O3 -ffast-math -fopenmp -c -o $@ $<
	gcc -DCFODE_STATIC -std=c99 -g -O0 -ffast-math -fopenmp -S -o $(<:.c=.s) $<



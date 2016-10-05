
#make for sample use of lsoda() stiff ode solver

#CFLAGS= -g -I$(HOME)/libs/
CFLAGS=  -g

LIB=-L$(HOME)/libs

FILES= daxpy.o ddot.o dgefa.o dgesl.o dscal.o idamax.o lsoda.o
#FILES=
#SAMPL= example.o
SAMPL= react1.o

example: $(SAMPL) $(FILES)
	cc -o example $(SAMPL) $(FILES) -s $(LIB) -lq -lm 

libq.a: $(FILES)
	ar cr libq.a `lorder $(FILES) | tsort`; ranlib libq.a


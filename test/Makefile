include ../Make.inc

CFLAGS-add += -I../src -DNDEBUG
LDFLAGS-add += -L../src -llsoda -lm

all: test

test.o: test.c
	$(CC) $(CPPFLAGS) $(CFLAGS-add) $(CFLAGS) $(fPIC) -c $< -o $@

test: test.o
	$(CC) $(LDFLAGS)  -o $@ $^ $(LDFLAGS-add)

clean:
	rm -f test.o test

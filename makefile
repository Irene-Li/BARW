# makefile for bws project
GIT_VERSION := $(shell git describe --abbrev=4 --always --tags)

CC = gcc
LDFLAGS = -lm 
CFLAGS = -Wall -DVERSION=\"$(GIT_VERSION)\"
OPTIM = -O3
BIN = bws

all: ${BIN}

clean:
	rm -f ${BIN} *.o

git:
	git add .
	git commit -m "$m"
	git push
	git push --tags
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) $(OPTIM) -o $(BIN)

bws_OBJS = bws.o rng/mt19937ar_clean_bkp.o lattice_core.o
bws: ${bws_OBJS}
	$(CC) $(CFLAGS) $(OPTIM) -o bws ${bws_OBJS} $(LDFLAGS)

.SUFFIXES : .o .c .h
.c.o :
	${CC} ${OPTIM} ${CFLAGS} -c $< -o $@


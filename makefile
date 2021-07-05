# makefile for bws project
GIT_VERSION := $(shell git describe --abbrev=4 --always --tags)

CC = gcc
LDFLAGS = -lm 
CFLAGS = -Wall -DVERSION=\"$(GIT_VERSION)\"
OPTIM = -O3
BIN = bws
SRC = bws.c 

all:
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) $(OPTIM) -o $(BIN)
git:
	git add .
	git commit -m "$m"
	git push
	git push --tags


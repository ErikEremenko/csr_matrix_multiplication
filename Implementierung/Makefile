### CONFIGURATION ###
WARNINGS := -Wall -Wextra
OPTIMIZATION := -O3
VERSION := -std=c17
AVX := -mavx
#####################

CC := gcc
CFLAGS := $(WARNINGS) $(OPTIMIZATION)

SOURCES = main.c matrixutils.c constants.c utils.c matrix.c

main: main.c
	$(CC) $(SOURCES) $(CFLAGS) $(VERSION) $(AVX) -o main

clean:
	rm -f main

CFLAGS = -O3 -g -lm -fopenmp
SRC = main.c util.c util.h constants.c constants.h tree_code.h
CC = gcc

glens: main.o util.o constants.o tree_code.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf ./*.o glens

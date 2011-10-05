CFLAGS = -O3 -g -lm
SRC = main.c tree_struct.c tree_struct.h util.c util.h constants.c constants.h
CC = gcc

glens: main.o util.o constants.o tree_struct.o
		$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf ./*.o glens

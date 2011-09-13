CFLAGS = -O3 -g -lm
SRC = main.c util.c util.h constants.c constants.h
CC = gcc

glens: main.o util.o constants.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf ./*.o glens

CFLAGS = -g -lm
SRC = constants.c constants.h lens_gen.c util.c util.h
CC = gcc

lens_gen: main.o constants.o util.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf ./*.o lens_gen

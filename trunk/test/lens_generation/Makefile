CFLAGS = -g -lm
SRC = main.c util.c util.h constants.c constants.h lens_gen.c global.h util.h
CC = gcc

lens_gen: main.o constants.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf ./*.o glens

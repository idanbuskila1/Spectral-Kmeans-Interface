CC=gcc -ansi -g
CFLAGS= -Wall -Wextra -Werror -pedantic-errors -lm
DEPS = spkmeans.h
OBJ = spkmeans.o jacobi.o kmeans.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

spkmeans: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
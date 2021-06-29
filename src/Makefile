CC=gcc

CFLAGS = -std=c11 -Wall -Wextra -Wshadow -Wpedantic -g -O0 -fsanitize=leak

LDLIBS= -lm

ALL = projet
all : $(ALL)

projet   : projet.c

fmt:
	clang-format -i -style=file projet.c
	clang-format -i -style=file projet.h

clean:
	rm -rf *.o
cleanexe:
	rm -rf $(ALL)

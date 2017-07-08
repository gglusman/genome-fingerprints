
#CFLAGS=-g
CFLAGS=-O3
#LFLAGS=-g
LFLAGS=-O3
CC=/bin/gcc

fpc:	fpc.o
	${CC} ${LFLAGS} -o fpc fpc.o -lm
	#cp fpc ~/bin/

fpc.o:	fpc.c


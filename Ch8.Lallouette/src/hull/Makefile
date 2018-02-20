# to install executable into $(BINDIR),
# and library into $(LIBDIR),
#	type "make".

#CC	= cc
CC	= gcc
AR	= ar	#/usr/ccs/bin/ar for Solaris
CFLAGS	= 
OBJS	= hull.o ch.o io.o rand.o pointops.o fg.o
HDRS	= hull.h points.h pointsites.h stormacs.h
SRC	= hull.c ch.c io.c rand.c pointops.c fg.c
BINDIR	= ../bin
LIBDIR	= .
LIB	= $(LIBDIR)/libhull.a

all	: $(OBJS)
	$(AR) -rs $(LIB) $(OBJS)


# the following probably needn't be changed

$(OBJS) : $(HDRS)
	$(CC) -c $(CFLAGS) $(SRC)

hullmain.o	: $(HDRS)

$(PROG)	: $(OBJS) hullmain.o
	$(CC) $(CFLAGS) $(OBJS) hullmain.o -o $(PROG) -lm

rsites$(EXE)	: rsites.c rand.o
	$(CC) $(CFLAGS) rand.o -o rsites rsites.c -lm

clean	:
	-rm -f $(OBJS) $(LIB) hullmain.o core a.out

# DO NOT DELETE

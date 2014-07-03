# Makefile for conversiontools

NAME	= conversiontools
VERSION	= $(shell git describe --tags --long)

CC		= gcc
CFLAGS	= -O3 -mcmodel=medium -Wall -pedantic -I$(LOCAL_LIB_PATH)/include -DVERSION=\"${VERSION}\"
LIBS	= -L$(LOCAL_LIB_PATH)/lib -lm -liof -lsilo -lartsfc -lgicreader

SRCS	= $(wildcard *.c)
EXES	= $(SRCS:.c=)

# Rules

all: $(EXES)

.c:
	$(CC) $(CFLAGS) -DNAME=\"$@\" $@.c $(LIBS) -o $@

clean:
	rm -f *~ *.o $(EXES)

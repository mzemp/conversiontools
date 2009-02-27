# Names

NAME	= conversiontools
TOOLS	= ts2ta ts2tb ts2gb ta2ts tb2ts gb2ts\
	tsspp2tsdpp tsdpp2tsspp tscom2tsphy ts2silo\
	aa2as as2aa cas eas gicb2ts artb2ts
EXT	= 64
VERSION = 1.1

# Compiler stuff

CC	= gcc
CFLAGS	= -O3 -Wall -I$(LOCAL_LIB_PATH)/include
LIBS	= -L$(LOCAL_LIB_PATH)/lib -lm -liof -lsilo -lgic_reader

# Rules

all:    $(TOOLS)

.c:
	$(CC) $(CFLAGS) $@.c $(LIBS) -o $@
ifneq ($(EXT),)
	mv $@ $@$(EXT)
endif

clean:
	-rm -f *.o *~
ifneq ($(EXT),)
	-rm *$(EXT)
else
	-rm $(TOOLS)
endif

tar:
	cd ..; tar cvf - $(NAME)/Makefile $(NAME)/*.c $(NAME)/*.h > $(NAME)-$(VERSION).tar


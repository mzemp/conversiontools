# Names

NAME	= conversiontools
TOOLS	= tipsy_xdr_2_tipsy_ascii \
	tipsy_xdr_2_tipsy_nb \
	tipsy_xdr_2_gadget_nb \
	tipsy_xdr_2_silo \
	tipsy_ascii_2_tipsy_xdr \
	tipsy_nb_2_tipsy_xdr \
	gadget_nb_2_tipsy_xdr \
	art_nb_2_tipsy_xdr \
	gic_nb_2_tipsy_xdr \
	tipsy_xdr_spp_2_tipsy_xdr_dpp \
	tipsy_xdr_dpp_2_tipsy_xdr_spp \
	tipsy_xdr_comoving_2_tipsy_xdr_physical \
	array_xdr_2_array_ascii \
	array_ascii_2_array_xdr \
	array_art_nb_2_array_xdr \
	combine_array_xdr \
	extract_array_xdr
EXT	= 
VERSION = 1.9

# Compiler stuff

CC	= gcc
CFLAGS	= -O3 -mcmodel=medium -Wall -I$(LOCAL_LIB_PATH)/include
LIBS	= -L$(LOCAL_LIB_PATH)/lib -lm -liof -lsilo -lart_sfc -lgic_reader

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


/* 
** tipsy_xdr_2_tipsy_xdr.c
**
** Written by Marcel Zemp
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <iof.h>

void usage(void);

int main(int argc, char **argv) {

	int i;
	int positionprecision, verboselevel;
	int writegas, writedark, writestar;
	TIPSY_HEADER thin, thout;
	TIPSY_GAS_PARTICLE tgp;
	TIPSY_DARK_PARTICLE tdp;
	TIPSY_STAR_PARTICLE tsp;
	TIPSY_GAS_PARTICLE_DPP tgpdpp;
	TIPSY_DARK_PARTICLE_DPP tdpdpp;
	TIPSY_STAR_PARTICLE_DPP tspdpp;
	XDR xdrsin, xdrsout;

	positionprecision = 0;
	verboselevel = 0;
	writegas = 1;
	writedark = 1;
	writestar = 1;
	i = 1;
	while (i < argc) {
		if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
			usage();
			}
		if (strcmp(argv[i],"-version") == 0) {
			fprintf(stderr,"%s (%s)\n",NAME,VERSION);
			exit(1);
			}
		if (strcmp(argv[i],"-spp") == 0) {
			positionprecision = 0;
			i++;
			}
		else if (strcmp(argv[i],"-dpp") == 0) {
			positionprecision = 1;
			i++;
			}
		else if (strcmp(argv[i],"-writegas") == 0) {
			i++;
			if (i >= argc) usage();
			writegas = atoi(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-writedark") == 0) {
			i++;
			if (i >= argc) usage();
			writedark = atoi(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-writestar") == 0) {
			i++;
			if (i >= argc) usage();
			writestar = atoi(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-verbose") == 0) {
			verboselevel = 1;
			i++;
			}
		else {
			usage();
			}
		}
	xdrstdio_create(&xdrsin,stdin,XDR_DECODE);
	read_tipsy_xdr_header(&xdrsin,&thin);
	/*
	** Recalculate header
	*/
	thout = thin;
	if (writegas == 0) thout.ngas = 0;
	if (writedark == 0) thout.ndark = 0;
	if (writestar == 0) thout.nstar = 0;
	thout.ntotal = thout.ngas+thout.ndark+thout.nstar;
	/*
	** Write new file
	*/
	xdrstdio_create(&xdrsout,stdout,XDR_ENCODE);
	write_tipsy_xdr_header(&xdrsout,&thout);
	if (positionprecision == 0) {
		for (i = 0; i < thin.ngas; i++) {
			read_tipsy_xdr_gas(&xdrsin,&tgp);
			if (writegas == 1) write_tipsy_xdr_gas(&xdrsout,&tgp);
			}
		for (i = 0; i < thin.ndark; i++) {
			read_tipsy_xdr_dark(&xdrsin,&tdp);
			if (writedark == 1) write_tipsy_xdr_dark(&xdrsout,&tdp);
			}
		for (i = 0; i < thin.nstar; i++) {
			read_tipsy_xdr_star(&xdrsin,&tsp);
			if (writestar == 1) write_tipsy_xdr_star(&xdrsout,&tsp);
			}
		}
	if (positionprecision == 1) {
		for (i = 0; i < thin.ngas; i++) {
			read_tipsy_xdr_gas_dpp(&xdrsin,&tgpdpp);
			if (writegas == 1) write_tipsy_xdr_gas_dpp(&xdrsout,&tgpdpp);
			}
		for (i = 0; i < thin.ndark; i++) {
			read_tipsy_xdr_dark_dpp(&xdrsin,&tdpdpp);
			if (writedark == 1) write_tipsy_xdr_dark_dpp(&xdrsout,&tdpdpp);
			}
		for (i = 0; i < thin.nstar; i++) {
			read_tipsy_xdr_star_dpp(&xdrsin,&tspdpp);
			if (writestar == 1) write_tipsy_xdr_star_dpp(&xdrsout,&tspdpp);
			}
		}
	xdr_destroy(&xdrsin);
	xdr_destroy(&xdrsout);
	if (verboselevel >= 0) {
		fprintf(stderr,"Time: %g Ntotal: %u Ngas: %u Ndark: %u Nstar: %u\n",
			thout.time,thout.ntotal,thout.ngas,thout.ndark,thout.nstar);
		}
	exit(0);
	}

void usage(void) {

	fprintf(stderr,"\n");
	fprintf(stderr,"%s (%s)\n",NAME,VERSION); 
	fprintf(stderr,"\n");
	fprintf(stderr,"Program converts tipsy XDR format to tipsy XDR format.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Please specify the following parameters:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-spp               : set this flag if input and output files have single precision positions (default)\n");
	fprintf(stderr,"-dpp               : set this flag if input and output files have double precision positions\n");
	fprintf(stderr,"-writegas <value>  : 0 = don't write out gas / 1 = write out gas (default: 1)\n");
	fprintf(stderr,"-writedark <value> : 0 = don't write out dark matter / 1 = write out dark matter (default: 1)\n");
	fprintf(stderr,"-writestar <value> : 0 = don't write out stars / 1 = write out stars (default: 1)\n");
	fprintf(stderr,"-verbose : verbose\n");
	fprintf(stderr,"< <name>           : input file in tipsy XDR format\n");
	fprintf(stderr,"> <name>           : output file in tipsy XDR format\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Other options:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-h or -help : display this help and exit\n");
	fprintf(stderr,"-version    : display version information and exit\n");
	fprintf(stderr,"\n");
	exit(1);
	}

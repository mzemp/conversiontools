/* 
** tipsy_nb_2_tipsy_xdr.c
**
** Written by Marcel Zemp
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <iof.h>

void usage(void);

int main(int argc, char **argv) {

	int i;
	int positionprecision, verboselevel;
	TIPSY_HEADER th;
	TIPSY_GAS_PARTICLE tgp;
	TIPSY_DARK_PARTICLE tdp;
	TIPSY_STAR_PARTICLE tsp;
	TIPSY_GAS_PARTICLE_DPP tgpdpp;
	TIPSY_DARK_PARTICLE_DPP tdpdpp;
	TIPSY_STAR_PARTICLE_DPP tspdpp;
	XDR xdrs;

	positionprecision = 0;
	verboselevel = 0;
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
		else if (strcmp(argv[i],"-verbose") == 0) {
			verboselevel = 1;
			i++;
			}
		else {
			usage();
			}
		}
	xdrstdio_create(&xdrs,stdout,XDR_DECODE);
	read_tipsy_nb_header(stdin,&th);
	write_tipsy_xdr_header(&xdrs,&th);
	if (positionprecision == 0) {
		for (i = 0; i < th.ngas; i++) {
			read_tipsy_nb_gas(stdin,&tgp);
			write_tipsy_xdr_gas(&xdrs,&tgp);
			}
		for (i = 0; i < th.ndark; i++) {
			read_tipsy_nb_dark(stdin,&tdp);
			write_tipsy_xdr_dark(&xdrs,&tdp);
			}
		for (i = 0; i < th.nstar; i++) {
			read_tipsy_nb_star(stdin,&tsp);
			write_tipsy_xdr_star(&xdrs,&tsp);
			}
		}
	if (positionprecision == 1) {
		for (i = 0; i < th.ngas; i++) {
			read_tipsy_nb_gas_dpp(stdin,&tgpdpp);
			write_tipsy_xdr_gas_dpp(&xdrs,&tgpdpp);
			}
		for (i = 0; i < th.ndark; i++) {
			read_tipsy_nb_dark_dpp(stdin,&tdpdpp);
			write_tipsy_xdr_dark_dpp(&xdrs,&tdpdpp);
			}
		for (i = 0; i < th.nstar; i++) {
			read_tipsy_nb_star_dpp(stdin,&tspdpp);
			write_tipsy_xdr_star_dpp(&xdrs,&tspdpp);
			}
		}
	xdr_destroy(&xdrs);
	if (verboselevel >= 0) {
		fprintf(stderr,"Time: %g Ntotal: %u Ngas: %u Ndark: %u Nstar: %u\n",
			th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
		}
	exit(0);
	}

void usage(void) {

	fprintf(stderr,"\n");
	fprintf(stderr,"%s (%s)\n",NAME,VERSION);
	fprintf(stderr,"\n");
	fprintf(stderr,"Program converts tipsy native binary format to tipsy XDR format.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Please specify the following parameters:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-spp     : set this flag if input and output files have single precision positions (default)\n");
	fprintf(stderr,"-dpp     : set this flag if input and output files have double precision positions\n");
	fprintf(stderr,"-verbose : verbose\n");
	fprintf(stderr,"< <name> : input file in tipsy native binary format\n");
	fprintf(stderr,"> <name> : output file in tipsy XDR format\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Other options:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-h or -help : display this help and exit\n");
	fprintf(stderr,"-version    : display version information and exit\n");
	fprintf(stderr,"\n");
	exit(1);
	}

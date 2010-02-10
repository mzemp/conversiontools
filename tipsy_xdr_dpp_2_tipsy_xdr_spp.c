/* 
** tsdpp2tsspp.c
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
    int verboselevel;
    TIPSY_HEADER th;
    TIPSY_GAS_PARTICLE tgp;
    TIPSY_DARK_PARTICLE tdp;
    TIPSY_STAR_PARTICLE tsp;
    TIPSY_GAS_PARTICLE_DPP tgpdpp;
    TIPSY_DARK_PARTICLE_DPP tdpdpp;
    TIPSY_STAR_PARTICLE_DPP tspdpp;
    XDR xdrsin, xdrsout;

    verboselevel = 0;
    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-v") == 0) {
	    verboselevel = 1;
	    i++;
	    }
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    xdrstdio_create(&xdrsin,stdin,XDR_DECODE);
    xdrstdio_create(&xdrsout,stdout,XDR_ENCODE);
    read_tipsy_xdr_header(&xdrsin,&th);
    write_tipsy_xdr_header(&xdrsout,&th);
    for (i = 0; i < th.ngas; i++) {
	read_tipsy_xdr_gas_dpp(&xdrsin,&tgpdpp);
	copy_tgpdpp_to_tgp(&tgpdpp,&tgp);
	write_tipsy_xdr_gas(&xdrsout,&tgp);
	}
    for (i = 0; i < th.ndark; i++) {
	read_tipsy_xdr_dark_dpp(&xdrsin,&tdpdpp);
	copy_tdpdpp_to_tdp(&tdpdpp,&tdp);
	write_tipsy_xdr_dark(&xdrsout,&tdp);
	}
    for (i = 0; i < th.nstar; i++) {
	read_tipsy_xdr_star_dpp(&xdrsin,&tspdpp);
	copy_tspdpp_to_tsp(&tspdpp,&tsp);
	write_tipsy_xdr_star(&xdrsout,&tsp);
	}
    xdr_destroy(&xdrsin);
    xdr_destroy(&xdrsout);
    if (verboselevel >= 0) {
	fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
		th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
	}
    exit(0);
    }

void usage(void) {
 
    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts tipsy XDR format with double precision positions\n");
    fprintf(stderr,"to tipsy XDR format with single precision positions.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"< <name> : input file in tipsy XDR format with double precision positions\n");
    fprintf(stderr,"> <name> : output file in tipsy XDR format with single precision positions\n");
    fprintf(stderr,"\n");
    exit(1);
    }

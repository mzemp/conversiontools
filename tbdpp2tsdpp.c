/* 
** tbdpp2tsdpp.c
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <iof.h>

void usage(void);

int main(int argc, char **argv) {

    int i;
    TIPSY_HEADER th;
    GAS_PARTICLE_DPP gpdpp;
    DARK_PARTICLE_DPP dpdpp;
    STAR_PARTICLE_DPP spdpp;
    XDR xdrs;

    i = 1;
    while (i < argc) {
	if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    xdrstdio_create(&xdrs,stdout,XDR_DECODE);
    read_tipsy_binary_header(stdin,&th);
    write_tipsy_standard_header(&xdrs,&th);
    for (i = 0; i < th.ngas; i++) {
	read_tipsy_binary_gas_dpp(stdin,&gpdpp);
	write_tipsy_standard_gas_dpp(&xdrs,&gpdpp);
	}
    for (i = 0; i < th.ndark; i++) {
	read_tipsy_binary_dark_dpp(stdin,&dpdpp);
	write_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
	}
    for (i = 0; i < th.nstar; i++) {
	read_tipsy_binary_star_dpp(stdin,&spdpp);
	write_tipsy_standard_star_dpp(&xdrs,&spdpp);
	}
    xdr_destroy(&xdrs);
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts tipsy binary format to tipsy standard binary format\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"< <name> : input file in tipsy binary format with double precision positions\n");
    fprintf(stderr,"> <name> : output file in tipsy standard binary format with double precision positions\n");
    fprintf(stderr,"\n");
    exit(1);
    }

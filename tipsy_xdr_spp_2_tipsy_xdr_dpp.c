/* 
** tsspp2tsdpp.c
**
** written by Marcel Zemp
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
    TIPSY_HEADER th;
    GAS_PARTICLE gp;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
    GAS_PARTICLE_DPP gpdpp;
    DARK_PARTICLE_DPP dpdpp;
    STAR_PARTICLE_DPP spdpp;
    XDR xdrsin, xdrsout;

    i = 1;
    while (i < argc) {
	if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    xdrstdio_create(&xdrsin,stdin,XDR_DECODE);
    xdrstdio_create(&xdrsout,stdout,XDR_ENCODE);
    read_tipsy_standard_header(&xdrsin,&th);
    write_tipsy_standard_header(&xdrsout,&th);
    for (i = 0; i < th.ngas; i++) {
	read_tipsy_standard_gas(&xdrsin,&gp);
	copy_gp_to_gpdpp(&gp,&gpdpp);
	write_tipsy_standard_gas_dpp(&xdrsout,&gpdpp);
	}
    for (i = 0; i < th.ndark; i++) {
	read_tipsy_standard_dark(&xdrsin,&dp);
	copy_dp_to_dpdpp(&dp,&dpdpp);
	write_tipsy_standard_dark_dpp(&xdrsout,&dpdpp);
	}
    for (i = 0; i < th.nstar; i++) {
	read_tipsy_standard_star(&xdrsin,&sp);
	copy_sp_to_spdpp(&sp,&spdpp);
	write_tipsy_standard_star_dpp(&xdrsout,&spdpp);
	}
    xdr_destroy(&xdrsin);
    xdr_destroy(&xdrsout);
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
    exit(0);
    }

void usage(void) {
 
    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts tipsy standard binary format with single precision positions\n");
    fprintf(stderr,"to tipsy standard binary format with double precision positions.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"< <name> : input file in tipsy standard binary format with single precision positions\n");
    fprintf(stderr,"> <name> : output file in tipsy standard binary format with double precision positions\n");
    fprintf(stderr,"\n");
    exit(1);
    }

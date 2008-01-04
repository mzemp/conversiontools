/* 
** tb2ts.c
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "IOfunctions.h"

void usage(void);

int main(int argc, char **argv) {

    int i;
    TIPSY_HEADER th;
    GAS_PARTICLE gp;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
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
    assert(fread(&th,sizeof(TIPSY_HEADER),1,stdin) == 1);
    write_tipsy_standard_header(&xdrs,&th);
    for (i = 0; i < th.ngas; i++) {
	assert(fread(&gp,sizeof(GAS_PARTICLE),1,stdin) == 1);
	write_tipsy_standard_gas(&xdrs,&gp);
	}
    for (i = 0; i < th.ndark; i++) {
	assert(fread(&dp,sizeof(DARK_PARTICLE),1,stdin) == 1);
	write_tipsy_standard_dark(&xdrs,&dp);
	}
    for (i = 0; i < th.nstar; i++) {
	assert(fread(&sp,sizeof(STAR_PARTICLE),1,stdin) == 1);
	write_tipsy_standard_star(&xdrs,&sp);
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
    fprintf(stderr,"Please specify the following parametes:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"< input file in tipsy binary format\n");
    fprintf(stderr,"> output file in tipsy standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

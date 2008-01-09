/* 
** tadpp2tsdpp.c
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
    TIPSY_STRUCTURE_DPP *tsdpp;

    i = 1;
    while (i < argc) {
	if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    tsdpp = malloc(sizeof(TIPSY_STRUCTURE));
    assert(tsdpp != NULL);
    tsdpp->th = NULL;
    tsdpp->gpdpp = NULL;
    tsdpp->dpdpp = NULL;
    tsdpp->spdpp = NULL;
    read_tipsy_ascii_dpp(stdin,tsdpp);
    write_tipsy_standard_dpp(stdout,tsdpp);
    fprintf(stderr,"\n");
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    tsdpp->th->time,tsdpp->th->ntotal,tsdpp->th->ngas,tsdpp->th->ndark,tsdpp->th->nstar);
    fprintf(stderr,"\n");
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts tipsy ascii format to tipsy standard binary format\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"< <name> : input file in tipsy ascii format with double precision positions\n");
    fprintf(stderr,"> <name> : output file in tipsy standard binary format with double precision positions\n");
    fprintf(stderr,"\n");
    exit(1);
    }

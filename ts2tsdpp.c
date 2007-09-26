/* 
** ts2tsdpp.c
**
** Program written in order to convert tipsy standard binary format with single precision positions
** to tipsy standard binary format with double precision positions
**
** written by Marcel Zemp, mzemp@ucolick.org, April 2007
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "IOfunctions.h"

void usage(void);

int main(int argc, char **argv) {

    int i;
    TIPSY_STRUCTURE *ts;

    i = 1;
    while (i < argc) {
	if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    ts = malloc(sizeof(TIPSY_STRUCTURE));
    assert(ts != NULL);
    ts->th = NULL;
    ts->gp = NULL;
    ts->dp = NULL;
    ts->sp = NULL;
    read_tipsy_standard(stdin,ts);
    write_tipsy_standard(stdout,ts);
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    ts->th->time,ts->th->ntotal,ts->th->ngas,ts->th->ndark,ts->th->nstar);
    exit(0);
    }

void usage(void) {
 
    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts tipsy standard binary format with single precision positions\n");
    fprintf(stderr,"to tipsy standard binary format with double precision positions\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parametes:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"< input file in tipsy standard binary format with single precision positions\n");
    fprintf(stderr,"> output file in tipsy standard binary format with double precision positions\n");
    fprintf(stderr,"\n");
    exit(1);
    }

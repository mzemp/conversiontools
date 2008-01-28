/* 
** ts2ta.c
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
    int positionprecision;
    TIPSY_HEADER *th;
    TIPSY_STRUCTURE *ts;
    TIPSY_STRUCTURE_DPP *tsdpp;

    th = NULL;
    ts = NULL;
    tsdpp = NULL;
    positionprecision = 0;
    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-spp") == 0) {
            positionprecision = 0;
            i++;
            }
	else if (strcmp(argv[i],"-dpp") == 0) {
            positionprecision = 1;
            i++;
            }
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    if (positionprecision == 0) {
	tsdpp = NULL;
	ts = malloc(sizeof(TIPSY_STRUCTURE));
	assert(ts != NULL);
	ts->th = NULL;
	ts->gp = NULL;
	ts->dp = NULL;
	ts->sp = NULL;
	read_tipsy_standard(stdin,ts);
	write_tipsy_ascii(stdout,ts);
	th = ts->th;
	}
    else if (positionprecision == 1) {
	ts = NULL;
	tsdpp = malloc(sizeof(TIPSY_STRUCTURE_DPP));
	assert(tsdpp != NULL);
	tsdpp->th = NULL;
	tsdpp->gpdpp = NULL;
	tsdpp->dpdpp = NULL;
	tsdpp->spdpp = NULL;
	read_tipsy_standard_dpp(stdin,tsdpp);
	write_tipsy_ascii_dpp(stdout,tsdpp);
	th = tsdpp->th;
	}
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    th->time,th->ntotal,th->ngas,th->ndark,th->nstar);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts tipsy standard binary format to tipsy ascii format\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp     : set this flag if input and output file have single precision positions (default)\n");
    fprintf(stderr,"-dpp     : set this flag if input and output file have double precision positions\n");
    fprintf(stderr,"< <name> : input file in tipsy standard binary format\n");
    fprintf(stderr,"> <name> : output file in tipsy ascii format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

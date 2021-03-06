/* 
** tipsy_xdr_2_tipsy_ascii.c
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
	TIPSY_HEADER *th;
	TIPSY_STRUCTURE *ts;
	TIPSY_STRUCTURE_DPP *tsdpp;

	th = NULL;
	ts = NULL;
	tsdpp = NULL;
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
	if (positionprecision == 0) {
		tsdpp = NULL;
		ts = malloc(sizeof(TIPSY_STRUCTURE));
		assert(ts != NULL);
		ts->th = NULL;
		ts->tgp = NULL;
		ts->tdp = NULL;
		ts->tsp = NULL;
		read_tipsy_xdr(stdin,ts);
		write_tipsy_ascii(stdout,ts);
		th = ts->th;
		}
	else if (positionprecision == 1) {
		ts = NULL;
		tsdpp = malloc(sizeof(TIPSY_STRUCTURE_DPP));
		assert(tsdpp != NULL);
		tsdpp->th = NULL;
		tsdpp->tgpdpp = NULL;
		tsdpp->tdpdpp = NULL;
		tsdpp->tspdpp = NULL;
		read_tipsy_xdr_dpp(stdin,tsdpp);
		write_tipsy_ascii_dpp(stdout,tsdpp);
		th = tsdpp->th;
		}
	if (verboselevel > 0) {
		fprintf(stderr,"Time: %g Ntotal: %u Ngas: %u Ndark: %u Nstar: %u\n",
			th->time,th->ntotal,th->ngas,th->ndark,th->nstar);
		}
	exit(0);
	}

void usage(void) {

	fprintf(stderr,"\n");
	fprintf(stderr,"%s (%s)\n",NAME,VERSION);
	fprintf(stderr,"\n");
	fprintf(stderr,"Program converts tipsy XDR format to tipsy ascii format.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Please specify the following parameters:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-spp     : set this flag if input and output files have single precision positions (default)\n");
	fprintf(stderr,"-dpp     : set this flag if input and output files have double precision positions\n");
	fprintf(stderr,"-verbose : verbose\n");
	fprintf(stderr,"< <name> : input file in tipsy XDR format\n");
	fprintf(stderr,"> <name> : output file in tipsy ascii format\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Other options:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-h or -help : display this help and exit\n");
	fprintf(stderr,"-version    : display version information and exit\n");
	fprintf(stderr,"\n");
	exit(1);
	}

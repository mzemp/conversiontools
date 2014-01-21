/* 
** array_xdr_2_array_ascii.c
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
	ARRAY_HEADER ah;
	ARRAY_PARTICLE ap;
	XDR xdrs;

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
	xdrstdio_create(&xdrs,stdin,XDR_DECODE);
	read_array_xdr_header(&xdrs,&ah);
	allocate_array_particle(&ah,&ap);
	assert(fprintf(stdout,"%d\n",ah.N[0]) > 0);
	if (ah.N[1] == 1) {
		assert(ah.N[2] == 0);
		assert(ah.N[3] == 0);
		for (i = 0; i < ah.N[0]; i++) {
			read_array_xdr_particle(&xdrs,&ah,&ap);
			assert(fprintf(stdout,"%d\n",ap.ia[0]) > 0);
			}
		}
	else if (ah.N[2] == 1) {
		assert(ah.N[1] == 0);
		assert(ah.N[3] == 0);
		for (i = 0; i < ah.N[0]; i++) {
			read_array_xdr_particle(&xdrs,&ah,&ap);
			assert(fprintf(stdout,"%.6e\n",ap.fa[0]) > 0);
			}
		}
	else if (ah.N[3] == 1) {
		assert(ah.N[1] == 0);
		assert(ah.N[2] == 0);
		for (i = 0; i < ah.N[0]; i++) {
			read_array_xdr_particle(&xdrs,&ah,&ap);
			assert(fprintf(stdout,"%.14e\n",ap.da[0]) > 0);
			}
		}
	xdr_destroy(&xdrs);
	if (verboselevel >= 0) {
		fprintf(stderr,"Ntotal: %d Ni: %d Nf: %d Nd: %d\n",
			ah.N[0],ah.N[1],ah.N[2],ah.N[3]);
		}
	exit(0);
	}

void usage(void) {

	fprintf(stderr,"\n");
	fprintf(stderr,"Program converts array in XDR format to array in ascii format.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Please specify the following parameters:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"< <name> : input file in XDR format\n");
	fprintf(stderr,"> <name> : output file in ascii format\n");
	fprintf(stderr,"\n");
	exit(1);
	}

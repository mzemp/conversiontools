/* 
** array_ascii_2_array_xdr.c
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

	for (i = 0; i < 4; i++) {
		ah.N[i] = 0;
		}
	verboselevel = 0;
	i = 1;
	while (i < argc) {
		if (strcmp(argv[i],"-i") == 0) {
			ah.N[1] = 1;
			ah.N[2] = 0;
			ah.N[3] = 0;
			i++;
			}
		else if (strcmp(argv[i],"-f") == 0) {
			ah.N[1] = 0;
			ah.N[2] = 1;
			ah.N[3] = 0;
			i++;
			}
		else if (strcmp(argv[i],"-d") == 0) {
			ah.N[1] = 0;
			ah.N[2] = 0;
			ah.N[3] = 1;
			i++;
			}
		else if (strcmp(argv[i],"-v") == 0) {
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
	if ((ah.N[1] == 0) && (ah.N[2] == 0) && (ah.N[3] == 0)) {
		fprintf(stderr,"No array type specified!\n");
		usage();
		}
	xdrstdio_create(&xdrs,stdout,XDR_ENCODE);
	allocate_array_particle(&ah,&ap);
	assert(fscanf(stdin,"%d",&ah.N[0]) == 1);
	write_array_xdr_header(&xdrs,&ah);
	if (ah.N[1] == 1) {
		for (i = 0; i < ah.N[0]; i++) {
			assert(fscanf(stdin,"%d",&ap.ia[0]) == 1);
			write_array_xdr_particle(&xdrs,&ah,&ap);
			}
		}
	else if (ah.N[2] == 1) {
		for (i = 0; i < ah.N[0]; i++) {
			assert(fscanf(stdin,"%f",&ap.fa[0]) == 1);
			write_array_xdr_particle(&xdrs,&ah,&ap);
			}
		}
	else if (ah.N[3] == 1) {
		for (i = 0; i < ah.N[0]; i++) {
			assert(fscanf(stdin,"%lf",&ap.da[0]) == 1);
			write_array_xdr_particle(&xdrs,&ah,&ap);
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
	fprintf(stderr,"Program converts array in ascii format to array in XDR format.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Please specify the following parameters:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-i		 : set this flag if it is an integer array\n");
	fprintf(stderr,"-f		 : set this flag if it is a float array\n");
	fprintf(stderr,"-d		 : set this flag if it is a double array\n");
	fprintf(stderr,"< <name> : input file in ascii format\n");
	fprintf(stderr,"> <name> : output file in XDR format\n");
	fprintf(stderr,"\n");
	exit(1);
	}

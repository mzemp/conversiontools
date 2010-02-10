/* 
** aartb2as.c
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
    int header, trailer;
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
    assert(fread(&header,sizeof(int),1,stdin) == 1);
    assert(fread(&ah.N[0],sizeof(int),1,stdin) == 1);
    assert(fread(&trailer,sizeof(int),1,stdin) == 1);
    assert(header == trailer);
    write_array_xdr_header(&xdrs,&ah);
    assert(fread(&header,sizeof(int),1,stdin) == 1);
    if (ah.N[1] == 1) {
	for (i = 0; i < ah.N[0]; i++) {
	    assert(fread(&ap.ia[0],sizeof(int),1,stdin) == 1);
	    write_array_xdr_particle(&xdrs,&ah,&ap);
	    }
	}
    else if (ah.N[2] == 1) {
	for (i = 0; i < ah.N[0]; i++) {
	    assert(fread(&ap.fa[0],sizeof(float),1,stdin) == 1);
	    write_array_xdr_particle(&xdrs,&ah,&ap);
	    }
	}
    else if (ah.N[3] == 1) {
	for (i = 0; i < ah.N[0]; i++) {
	    assert(fread(&ap.da[0],sizeof(double),1,stdin) == 1);
	    write_array_xdr_particle(&xdrs,&ah,&ap);
	    }
	}
    assert(fread(&trailer,sizeof(int),1,stdin) == 1);
    assert(header == trailer);
    xdr_destroy(&xdrs);
    if (verboselevel >= 0) {
	fprintf(stderr,"Ntotal: %d Ni: %d Nf: %d Nd: %d\n",
		ah.N[0],ah.N[1],ah.N[2],ah.N[3]);
	}
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts array in ART native binary format to array in XDR format.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-i       : set this flag if it is an integer array\n");
    fprintf(stderr,"-f       : set this flag if it is an float array\n");
    fprintf(stderr,"-d       : set this flag if it is an double array\n");
    fprintf(stderr,"< <name> : input file in ART native binary format\n");
    fprintf(stderr,"> <name> : output file in XDR format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

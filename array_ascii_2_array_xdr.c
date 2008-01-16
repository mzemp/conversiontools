/* 
** aa2as.c
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
    ARRAY_HEADER ah;
    ARRAY_PARTICLE ap;
    XDR xdrs;

    xdrstdio_create(&xdrs,stdout,XDR_ENCODE);
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
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    allocate_array_particle(&ah,&ap);
    assert(fscanf(stdin,"%d",&ah.N[0]) == 1);
    write_array_header(&xdrs,&ah);
    if (ah.N[1] == 1) {
	for (i = 0; i < ah.N[0]; i++) {
	    assert(fscanf(stdin,"%d",&ap.ia[0]) == 1);
	    write_array_particle(&xdrs,&ah,&ap);
	    }
	}
    else if (ah.N[2] == 1) {
	for (i = 0; i < ah.N[0]; i++) {
	    assert(fscanf(stdin,"%f",&ap.fa[0]) == 1);
	    write_array_particle(&xdrs,&ah,&ap);
	    }
	}
    else if (ah.N[3] == 1) {
	for (i = 0; i < ah.N[0]; i++) {
	    assert(fscanf(stdin,"%lf",&ap.da[0]) == 1);
	    write_array_particle(&xdrs,&ah,&ap);
	    }
	}
    xdr_destroy(&xdrs);
    fprintf(stderr,"Ntotal: %d Ni: %d Nf: %d Nd: %d\n",
	    ah.N[0],ah.N[1],ah.N[2],ah.N[3]);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts array in ascii format to array in standard binary format\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-i       : set this flag if it is an integer array\n");
    fprintf(stderr,"-f       : set this flag if it is an float array\n");
    fprintf(stderr,"-d       : set this flag if it is an double array\n");
    fprintf(stderr,"< <name> : input file in ascii format\n");
    fprintf(stderr,"> <name> : output file in standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

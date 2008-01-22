/* 
** eas.c
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
    int iindex = 0;
    int findex = 0;
    int dindex = 0;
    ARRAY_HEADER ah1, ah2;
    ARRAY_PARTICLE ap1, ap2;
    XDR xdrs1, xdrs2;

    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-i") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    iindex = atoi(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-f") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    findex = atoi(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-d") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    dindex = atoi(argv[i]);
	    i++;
	    }
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    xdrstdio_create(&xdrs1,stdin,XDR_DECODE);
    xdrstdio_create(&xdrs2,stdout,XDR_ENCODE);
    read_array_header(&xdrs1,&ah1);
    ah2.N[0] = ah1.N[0];
    if (iindex != 0) {
	assert(findex == 0);
	assert(dindex == 0);
	ah2.N[1] = 1;
	ah2.N[2] = 0;
	ah2.N[3] = 0;
	}
    if (findex != 0) {
	assert(iindex == 0);
	assert(dindex == 0);
	ah2.N[1] = 0;
	ah2.N[2] = 1;
	ah2.N[3] = 0;
	}
    if (dindex != 0) {
	assert(iindex == 0);
	assert(findex == 0);
	ah2.N[1] = 0;
	ah2.N[2] = 0;
	ah2.N[3] = 1;
	}
    write_array_header(&xdrs2,&ah2);
    allocate_array_particle(&ah1,&ap1);
    allocate_array_particle(&ah2,&ap2);

    for (i = 0; i < ah2.N[0]; i++) {
	read_array_particle(&xdrs1,&ah1,&ap1);
	if (iindex != 0) {
	    ap2.ia[0] = ap1.ia[iindex-1];
	    }
	else if (findex != 0) {
	    ap2.fa[0] = ap1.fa[findex-1];
	    }
	else if (dindex != 0) {
	    ap2.da[0] = ap1.da[dindex-1];
	    }
	write_array_particle(&xdrs2,&ah2,&ap2);
	}
    xdr_destroy(&xdrs1);
    xdr_destroy(&xdrs2);
    fprintf(stderr,"Ntotal_in: %d Ni_in: %d Nf_in: %d Nd_in: %d\n",
	    ah1.N[0],ah1.N[1],ah1.N[2],ah1.N[3]);
    fprintf(stderr,"Ntotal_out: %d Ni_out: %d Nf_out: %d Nd_out: %d\n",
	    ah2.N[0],ah2.N[1],ah2.N[2],ah2.N[3]);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program extracts specific array from array file in standard binary format\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-i <value> : index of integer array to be extracted or\n");
    fprintf(stderr,"-f <value> : index of float array to be extracted or\n");
    fprintf(stderr,"-d <value> : index of double array to be extracted\n");
    fprintf(stderr,"< <name>   : input file in standard binary format\n");
    fprintf(stderr,"> <name>   : output file in standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }
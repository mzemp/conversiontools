/* 
** cas.c
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

    int i, j;
    ARRAY_HEADER ah1, ah2, ah3;
    ARRAY_PARTICLE ap1, ap2, ap3;
    FILE *fp1, *fp2;
    XDR xdrs1, xdrs2, xdrs3;

    fp1 = NULL;
    fp2 = NULL;
    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-f1") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    fp1 = fopen(argv[i],"r");
	    i++;
	    }
	else if (strcmp(argv[i],"-f2") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    fp2 = fopen(argv[i],"r");
	    i++;
	    }
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
            usage();
            }
        else {
            usage();
            }
        }
    assert(fp1 != NULL);
    assert(fp2 != NULL);
    xdrstdio_create(&xdrs1,fp1,XDR_DECODE);
    xdrstdio_create(&xdrs2,fp2,XDR_DECODE);
    xdrstdio_create(&xdrs3,stdout,XDR_ENCODE);
    read_array_header(&xdrs1,&ah1);
    read_array_header(&xdrs2,&ah2);
    assert(ah1.N[0] == ah2.N[0]);
    ah3.N[0] = ah1.N[0];
    allocate_array_particle(&ah1,&ap1);
    allocate_array_particle(&ah2,&ap2);
    for (i = 1; i < 4; i++) {
	ah3.N[i] = ah1.N[i] + ah2.N[i];
	}
    allocate_array_particle(&ah3,&ap3);
    write_array_header(&xdrs3,&ah3);
    for (i = 0; i < ah3.N[0]; i++) {
	read_array_particle(&xdrs1,&ah1,&ap1);
	read_array_particle(&xdrs2,&ah2,&ap2);
	for (j = 0; j < ah1.N[1]; j++) {
	    ap3.ia[j] = ap1.ia[j];
	    }
	for (j = 0; j < ah2.N[1]; j++) {
	    ap3.ia[ah1.N[1]+j] = ap2.ia[j];
	    }
	for (j = 0; j < ah1.N[2]; j++) {
	    ap3.fa[j] = ap1.fa[j];
	    }
	for (j = 0; j < ah2.N[2]; j++) {
	    ap3.fa[ah1.N[2]+j] = ap2.fa[j];
	    }
	for (j = 0; j < ah1.N[3]; j++) {
	    ap3.da[j] = ap1.da[j];
	    }
	for (j = 0; j < ah2.N[3]; j++) {
	    ap3.da[ah1.N[3]+j] = ap2.da[j];
	    }
	write_array_particle(&xdrs3,&ah3,&ap3);
	}
    xdr_destroy(&xdrs1);
    xdr_destroy(&xdrs2);
    xdr_destroy(&xdrs3);
    fprintf(stderr,"Ntotal_1: %d Ni_1: %d Nf_1: %d Nd_1: %d\n",
	    ah1.N[0],ah1.N[1],ah1.N[2],ah1.N[3]);
    fprintf(stderr,"Ntotal_2: %d Ni_2: %d Nf_2: %d Nd_2: %d\n",
	    ah2.N[0],ah2.N[1],ah2.N[2],ah2.N[3]);
    fprintf(stderr,"Ntotal_out: %d Ni_out: %d Nf_out: %d Nd_out: %d\n",
	    ah3.N[0],ah3.N[1],ah3.N[2],ah3.N[3]);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program combines two array files in standard binary format.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-f1 <name> : input file 1 in standard binary format\n");
    fprintf(stderr,"-f2 <name> : input file 2 in standard binary format\n");
    fprintf(stderr,"> <name>   : output file in standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

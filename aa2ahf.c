/* 
** aa2ahf.c
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

    int i, Ntot;
    int record;
    int verboselevel;
    float value;

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
    assert(fscanf(stdin,"%d",&Ntot) == 1);
    record = sizeof(int);
    assert(fwrite(&record,sizeof(int),1,stdout) == 1);
    assert(fwrite(&Ntot,sizeof(int),1,stdout) == 1);
    assert(fwrite(&record,sizeof(int),1,stdout) == 1);
    record = Ntot*sizeof(float);
    assert(fwrite(&record,sizeof(int),1,stdout) == 1);
    for (i = 0; i < Ntot; i++) {
	assert(fscanf(stdin,"%f",&value) == 1);
	assert(fwrite(&value,sizeof(float),1,stdout) == 1);
	}
    assert(fwrite(&record,sizeof(int),1,stdout) == 1);
    if (verboselevel > 0) {
	fprintf(stderr,"Ntotal: %d\n",Ntot);
	}
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts array in ascii format to array in native binary format that can be read by hfind (Fortran style).\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"< <name> : input file in ascii format\n");
    fprintf(stderr,"> <name> : output file in native binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

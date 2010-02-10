/* 
** gadget_nb_2_tipsy_xdr.c
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
    double a, dx, dy, dz, dof, mmw, uvf;
    TIPSY_STRUCTURE *ts;

    a = 1;
    dx = 0;
    dy = 0;
    dz = 0;
    dof = 3;
    mmw = 1;
    uvf = 977.79219;
    verboselevel = 0;
    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-a") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    if (strcmp(argv[i],"time") == 0) {
		a = -1;
		}
	    else {
		a = atof(argv[i]);
		}
	    i++;
            }
	else if (strcmp(argv[i],"-drx") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    dx = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-dry") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    dy = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-drz") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    dz = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-dof") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    dof = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-mmw") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    mmw = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-uvf") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    uvf = atof(argv[i]);
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
    ts = malloc(sizeof(TIPSY_STRUCTURE));
    assert(ts != NULL);
    ts->th = NULL;
    ts->tgp = NULL;
    ts->tdp = NULL;
    ts->tsp = NULL;
    read_gadget_nb(stdin,ts,a,dx,dy,dz,dof,mmw,uvf);
    write_tipsy_xdr(stdout,ts);
    if (a == -1) {
	a = ts->th->time;
	}
    if (verboselevel >= 1) {
	fprintf(stderr,"Used values:\n");
	fprintf(stderr,"a   : %.6e\n",a);
	fprintf(stderr,"drx : %.6e LU\n",dx);
	fprintf(stderr,"dry : %.6e LU\n",dy);
	fprintf(stderr,"drz : %.6e LU\n",dz);
	fprintf(stderr,"dof : %.6e\n",dof);
	fprintf(stderr,"mmw : %.6e mp\n",mmw);
	fprintf(stderr,"uvf : %.6e m s^-1\n",uvf);
	}
    if (verboselevel >= 0) {
	fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
		ts->th->time,ts->th->ntotal,ts->th->ngas,ts->th->ndark,ts->th->nstar);
	}
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts gadget native binary format to tipsy XDR format.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-a <value>   : expansion factor [write -a time for cosmological runs] (default: 1)\n");
    fprintf(stderr,"-drx <value> : shift along x-axis [LU] (default: 0 LU)\n");
    fprintf(stderr,"-dry <value> : shift along y-axis [LU] (default: 0 LU)\n");
    fprintf(stderr,"-drz <value> : shift along z-axis [LU] (default: 0 LU)\n");
    fprintf(stderr,"-dof <value> : degrees of freedom of gas (default: 3 => gamma = (dof+2)/dof = 5/3)\n");
    fprintf(stderr,"-mmw <value> : mean molecular weight of gas [mp] (default: 1 mp)\n");
    fprintf(stderr,"-uvf <value> : internal unit of velocity [m s^-1] (default: 977.79219 m s^-1 => 977.79219 m s^-1 = 1 kpc Gyr^-1)\n");
    fprintf(stderr,"-v           : more informative output to screen\n");
    fprintf(stderr,"< <name>     : input file in gadget native binary format\n");
    fprintf(stderr,"> <name>     : output file in tipsy XDR format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

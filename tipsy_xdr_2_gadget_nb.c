/* 
** ts2gb.c
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "IOfunctions.h"

void usage(void);

int main(int argc, char **argv) {

    int i;
    double a, dx, dy, dz, dof, mmw, uvf;
    TIPSY_STRUCTURE *ts;

    a = 1;
    dx = 0;
    dy = 0;
    dz = 0;
    dof = 3;
    mmw = 1;
    uvf = 977.79219;

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
	else if (strcmp(argv[i],"-dx") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    dx = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-dy") == 0) {
	    i++;
	    if (i >= argc) {
                usage();
                }
	    dy = atof(argv[i]);
	    i++;
            }
	else if (strcmp(argv[i],"-dz") == 0) {
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
    ts->gp = NULL;
    ts->dp = NULL;
    ts->sp = NULL;
    read_tipsy_standard(stdin,ts);
    write_gadget_binary(stdout,ts,a,dx,dy,dz,dof,mmw,uvf);
    if (a == -1) {
	a = ts->th->time;
	}
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    ts->th->time,ts->th->ntotal,ts->th->ngas,ts->th->ndark,ts->th->nstar);
    fprintf(stderr,"Used values:\n");
    fprintf(stderr,"a   = %.6e\n",a);
    fprintf(stderr,"dx  = %.6e LU\n",dx);
    fprintf(stderr,"dy  = %.6e LU\n",dy);
    fprintf(stderr,"dz  = %.6e LU\n",dy);
    fprintf(stderr,"dof = %.6e\n",dof);
    fprintf(stderr,"mmw = %.6e mp\n",mmw);
    fprintf(stderr,"uvf = %.6e m s^-1\n",uvf);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts tipsy standard binary format to gadget binary format\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-a <value>   : expansion factor [write -a time for cosmological runs] (default: a = 1)\n");
    fprintf(stderr,"-dx <value>  : shift along x-axis [LU] (default: dx = 0 LU)\n");
    fprintf(stderr,"-dy <value>  : shift along y-axis [LU] (default: dy = 0 LU)\n");
    fprintf(stderr,"-dz <value>  : shift along z-axis [LU] (default: dz = 0 LU)\n");
    fprintf(stderr,"-dof <value> : degrees of freedom of gas (default: dof = 3 => gamma = (dof+2)/dof = 5/3)\n");
    fprintf(stderr,"-mmw <value> : mean molecular weight of gas [mp] (default: mmw = 1 mp)\n");
    fprintf(stderr,"-uvf <value> : internal unit of velocity [m s^-1] (default: uvf = 977.79219 m s^-1 => 977.79219 m s^-1 = 1 kpc Gyr^-1)\n");
    fprintf(stderr,"< <name>     : input file in tipsy standard binary format\n");
    fprintf(stderr,"> <name>     : output file in gadget binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

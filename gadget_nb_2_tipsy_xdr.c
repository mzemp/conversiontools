/* 
** gb2ts.c
**
** Program written in order to convert gadget binary format to tipsy standard binary format
**
** written by Marcel Zemp, mzemp@ucolick.org, May 2007
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
    read_gadget_binary(stdin,ts,a,dx,dy,dz,dof,mmw,uvf);
    write_tipsy_standard(stdout,ts);
    if (a == -1) {
	a = ts->th->time;
	}
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    ts->th->time,ts->th->ntotal,ts->th->ngas,ts->th->ndark,ts->th->nstar);
    fprintf(stderr,"a: %g dx: %g dy: %g dz: %g dof: %g mmw: %g uvf: %g\n",
	    a,dx,dy,dz,dof,mmw,uvf);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts gadget binary format to tipsy standard binary format\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parametes:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"< input file in gadget binary format\n");
    fprintf(stderr,"> output file in tipsy standard binary format\n");
    fprintf(stderr,"-a   : expansion factor [write -a time for cosmological rus] (default: a = 1)\n");
    fprintf(stderr,"-dx  : shift along x-axis (default: dx = 0)\n");
    fprintf(stderr,"-dy  : shift along y-axis (default: dy = 0)\n");
    fprintf(stderr,"-dz  : shift along z-axis (default: dz = 0)\n");
    fprintf(stderr,"-dof : degrees of freedom of gas (default: dof = 3 => gamma = (dof+2)/dof = 5/3)\n");
    fprintf(stderr,"-mmw : mean molecular weight of gas in units of proton mass (default: mmw = 1)\n");
    fprintf(stderr,"-uvf : internal unit of velocity in m s^-1 (default: uvf = 977.79219 => 977.79219 m s^-1 = 1 kpc Gyr^-1)\n");
    fprintf(stderr,"\n");
    exit(1);
    }

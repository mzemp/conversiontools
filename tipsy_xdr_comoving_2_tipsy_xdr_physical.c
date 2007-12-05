/* 
** comoving2physical.c
**
** Program written in order to convert output from cosmological simulations
** from comoving to physical coordinates
**
** written by Marcel Zemp, mzemp@ucolick.org, November 2007
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "IOfunctions.h"

void usage(void);

int main(int argc, char **argv) {

    int i, j;
    double rcen[3], vcen[3];
    double mscale, rscale, vscale, hubble;
    TIPSY_HEADER th;
    GAS_PARTICLE gp;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
    XDR xdrin, xdrout;

    for (j = 0; j < 3; j++) {
	rcen[j] = 0;
	vcen[j] = 0;
	}
    mscale = 1;
    rscale = 1;
    vscale = 1;
    hubble = 0;

    /*
    ** Read in arguments
    */

    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-rxcen") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    rcen[0] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rycen") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    rcen[1] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rzcen") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    rcen[2] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vxcen") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    vcen[0] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vycen") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    vcen[1] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vzcen") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    vcen[2] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-mscale") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    mscale = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rscale") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    rscale = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vscale") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    vscale = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Hubble") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    if (strcmp(argv[i],"sqrt_8pi_3") == 0) {
		hubble = sqrt(8*M_PI/3.0);
		}
	    else {
		hubble = atof(argv[i]);
		}
	    i++;
	    }
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
	    usage();
	    }
	else {
	    usage();
	    }
	}

    /*
    ** Read in old header and write out new one
    */

    xdrstdio_create(&xdrin,stdin,XDR_DECODE);
    read_tipsy_standard_header(&xdrin,&th);
    xdrstdio_create(&xdrout,stdout,XDR_ENCODE);
    write_tipsy_standard_header(&xdrout,&th);

    /*
    ** Add shifts to particles and write them out
    */

    for(i = 0; i < th.ngas; i++) {
	read_tipsy_standard_gas(&xdrin,&gp);
	for(j = 0; j < 3; j++) {
	    gp.pos[j] = (gp.pos[j]-rcen[j])*rscale;
	    gp.vel[j] = (gp.vel[j]-vcen[j])*vscale;
	    gp.vel[j] += hubble*gp.pos[j];
	    }
	gp.mass *= mscale;
	write_tipsy_standard_gas(&xdrout,&gp);
	}
    for(i = 0; i < th.ndark; i++) {
	read_tipsy_standard_dark(&xdrin,&dp);
	for(j = 0; j < 3; j++) {
	    dp.pos[j] = (dp.pos[j]-rcen[j])*rscale;
	    dp.vel[j] = (dp.vel[j]-vcen[j])*vscale;
	    dp.vel[j] += hubble*dp.pos[j];
	    }
	dp.mass *= mscale;
	write_tipsy_standard_dark(&xdrout,&dp);
	}
    for(i = 0; i < th.nstar; i++) {
	read_tipsy_standard_star(&xdrin,&sp);
	for(j = 0; j < 3; j++) {
	    sp.pos[j] = (sp.pos[j]-rcen[j])*rscale;
	    sp.vel[j] = (sp.vel[j]-vcen[j])*vscale;
	    sp.vel[j] += hubble*sp.pos[j];
	    }
	sp.mass *= mscale;
	write_tipsy_standard_star(&xdrout,&sp);
	}

    /*
    ** Clean up and write some output
    */

    xdr_destroy(&xdrin);
    xdr_destroy(&xdrout);

    fprintf(stderr,"\n");
    fprintf(stderr,"Command line:\n\n");
    for (i = 0; i < argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n\n");
    fprintf(stderr,"Used values:\n\n");
    fprintf(stderr,"rcen   = (%+.6e,%+.6e,%+.6e) LUold\n",rcen[0],rcen[1],rcen[2]);
    fprintf(stderr,"vcen   = (%+.6e,%+.6e,%+.6e) VUold\n",vcen[0],vcen[1],vcen[2]);
    fprintf(stderr,"mscale = %.6e MUnew/MUold\n",mscale);
    fprintf(stderr,"xscale = %.6e LUnew/LUold\n",rscale);
    fprintf(stderr,"vscale = %.6e VUnew/VUold\n",vscale);
    fprintf(stderr,"Hubble = %.6e VUnew/LUnew\n",hubble);
    fprintf(stderr,"\n");

    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Comoving2physical makes a coordinate transformation from comoving to physical coordinates.\n\n");
    fprintf(stderr,"m_new = m_old*mscale\n");
    fprintf(stderr,"r_new = (r_old-r_cen)*rscale\n");
    fprintf(stderr,"v_new = (v_old-v_cen)*vscale + Hubble*r_new\n\n");
    fprintf(stderr,"You can specify the following arguments:\n\n");
    fprintf(stderr,"-rxcen <value>  : x-coordinate of centre [LUold] (default: 0 LUold)\n");
    fprintf(stderr,"-rycen <value>  : y-coordinate of centre [LUold] (default: 0 LUold)\n");
    fprintf(stderr,"-rzcen <value>  : z-coordinate of centre [LUold] (default: 0 LUold)\n");
    fprintf(stderr,"-vxcen <value>  : x-velocity of centre [VUold] (default: 0 VUold)\n");
    fprintf(stderr,"-vycen <value>  : y-velocity of centre [VUold] (default: 0 VUold)\n");
    fprintf(stderr,"-vzcen <value>  : z-velocity of centre [VUold] (default: 0 VUold)\n");
    fprintf(stderr,"-mscale <value> : mass scale factor [MUnew/MUold] (default: 1 MUnew/MUold)\n");
    fprintf(stderr,"-rscale <value> : coordinate scale factor [LUnew/LUold] (default: 1 LUnew/LUold)\n");
    fprintf(stderr,"-vscale <value> : velocity scale factor [VUnew/VUold] (default: 1 VUnew/VUold)\n");
    fprintf(stderr,"-Hubble <value> : Hubble parameter [VUnew/LUnew] (default: 0 VUnew/LUnew; special value: sqrt_8pi_3)\n");
    fprintf(stderr,"< <name>        : name of input file (Tipsy standard binary)\n");
    fprintf(stderr,"> <name>        : name of output file (Tipsy standard binary)\n");
    fprintf(stderr,"\n");
    exit(1);
    }

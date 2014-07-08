/* 
** tipsy_xdr_2_gadget_nb.c
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
	uvf = ConversionFactors.kpc_per_Gyr_2_km_per_s*1e3;
	verboselevel = 0;
	i = 1;
	while (i < argc) {
		if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
			usage();
			}
		if (strcmp(argv[i],"-version") == 0) {
			fprintf(stderr,"%s (%s)\n",NAME,VERSION);
			exit(1);
			}
		if (strcmp(argv[i],"-a") == 0) {
			i++;
			if (i >= argc) usage();
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
			if (i >= argc) usage();
			dx = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-dry") == 0) {
			i++;
			if (i >= argc) usage();
			dy = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-drz") == 0) {
			i++;
			if (i >= argc) usage();
			dz = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-dof") == 0) {
			i++;
			if (i >= argc) usage();
			dof = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-mmw") == 0) {
			i++;
			if (i >= argc) usage();
			mmw = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-uvf") == 0) {
			i++;
			if (i >= argc) usage();
			uvf = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-verbose") == 0) {
			verboselevel = 1;
			i++;
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
	read_tipsy_xdr(stdin,ts);
	write_gadget_nb(stdout,ts,a,dx,dy,dz,dof,mmw,uvf);
	if (a == -1) {
		a = ts->th->time;
		}
	if (verboselevel > 0) {
		fprintf(stderr,"Used values:\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"a   : %.6e\n",a);
		fprintf(stderr,"drx : %.6e LU\n",dx);
		fprintf(stderr,"dry : %.6e LU\n",dy);
		fprintf(stderr,"drz : %.6e LU\n",dz);
		fprintf(stderr,"dof : %.6e\n",dof);
		fprintf(stderr,"mmw : %.6e mp\n",mmw);
		fprintf(stderr,"uvf : %.6e m s^-1\n",uvf);
		fprintf(stderr,"\n");
		fprintf(stderr,"Time: %g Ntotal: %u Ngas: %u Ndark: %u Nstar: %u\n",
			ts->th->time,ts->th->ntotal,ts->th->ngas,ts->th->ndark,ts->th->nstar);
		}
	exit(0);
	}

void usage(void) {

	double uvf = ConversionFactors.kpc_per_Gyr_2_km_per_s*1e3;

	fprintf(stderr,"\n");
	fprintf(stderr,"%s (%s)\n",NAME,VERSION);
	fprintf(stderr,"\n");
	fprintf(stderr,"Program converts tipsy XDR format to gadget native binary format.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Please specify the following parameters:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-a <value>   : expansion factor [write -a time for cosmological runs] (default: 1)\n");
	fprintf(stderr,"-drx <value> : shift along x-axis [LU] (default: 0 LU)\n");
	fprintf(stderr,"-dry <value> : shift along y-axis [LU] (default: 0 LU)\n");
	fprintf(stderr,"-drz <value> : shift along z-axis [LU] (default: 0 LU)\n");
	fprintf(stderr,"-dof <value> : degrees of freedom of gas (default: 3 => gamma = (dof+2)/dof = 5/3)\n");
	fprintf(stderr,"-mmw <value> : mean molecular weight of gas [mp] (default: 1 mp)\n");
	fprintf(stderr,"-uvf <value> : internal unit of velocity [m s^-1] (default: %g m s^-1 => %g m s^-1 = 1 kpc Gyr^-1)\n",uvf,uvf);
	fprintf(stderr,"-verbose     : verbose\n");
	fprintf(stderr,"< <name>     : input file in tipsy XDR format\n");
	fprintf(stderr,"> <name>     : output file in gadget native binary format\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Other options:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-h or -help : display this help and exit\n");
	fprintf(stderr,"-version    : display version information and exit\n");
	fprintf(stderr,"\n");
	exit(1);
	}

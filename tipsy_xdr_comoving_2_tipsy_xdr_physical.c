/* 
** tipsy_xdr_comoving_2_tipsy_xdr_physical.c
**
** Written by Marcel Zemp
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <iof.h>

void usage(void);

int main(int argc, char **argv) {

    int i, j;
    int positionprecision, verboselevel;
    double rcen[3], vcen[3];
    double mscale, rscale, vscale, hubble;
    TIPSY_HEADER th;
    TIPSY_GAS_PARTICLE tgp;
    TIPSY_DARK_PARTICLE tdp;
    TIPSY_STAR_PARTICLE tsp;
    TIPSY_GAS_PARTICLE_DPP tgpdpp;
    TIPSY_DARK_PARTICLE_DPP tdpdpp;
    TIPSY_STAR_PARTICLE_DPP tspdpp;
    XDR xdrsin, xdrsout;

    for (j = 0; j < 3; j++) {
	rcen[j] = 0;
	vcen[j] = 0;
	}
    mscale = 1;
    rscale = 1;
    vscale = 1;
    hubble = 0;
    positionprecision = 0;
    verboselevel = 0;
    /*
    ** Read in arguments
    */
    i = 1;
    while (i < argc) {
        if (strcmp(argv[i],"-spp") == 0) {
            positionprecision = 0;
            i++;
            }
        else if (strcmp(argv[i],"-dpp") == 0) {
            positionprecision = 1;
            i++;
            }
	else if (strcmp(argv[i],"-rxcen") == 0) {
	    i++;
	    if (i >= argc) usage();
	    rcen[0] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rycen") == 0) {
	    i++;
	    if (i >= argc) usage();
	    rcen[1] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rzcen") == 0) {
	    i++;
	    if (i >= argc) usage();
	    rcen[2] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vxcen") == 0) {
	    i++;
	    if (i >= argc) usage();
	    vcen[0] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vycen") == 0) {
	    i++;
	    if (i >= argc) usage();
	    vcen[1] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vzcen") == 0) {
	    i++;
	    if (i >= argc) usage();
	    vcen[2] = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-mscale") == 0) {
	    i++;
	    if (i >= argc) usage();
	    mscale = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rscale") == 0) {
	    i++;
	    if (i >= argc) usage();
	    rscale = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-vscale") == 0) {
	    i++;
	    if (i >= argc) usage();
	    vscale = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Hubble") == 0) {
	    i++;
	    if (i >= argc) usage();
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
    xdrstdio_create(&xdrsin,stdin,XDR_DECODE);
    xdrstdio_create(&xdrsout,stdout,XDR_ENCODE);
    read_tipsy_xdr_header(&xdrsin,&th);
    write_tipsy_xdr_header(&xdrsout,&th);
    /*
    ** Add shifts to particles and write them out
    */
    if (positionprecision == 0) {
	for(i = 0; i < th.ngas; i++) {
	    read_tipsy_xdr_gas(&xdrsin,&tgp);
	    for(j = 0; j < 3; j++) {
		tgp.pos[j] = (tgp.pos[j]-rcen[j])*rscale;
		tgp.vel[j] = (tgp.vel[j]-vcen[j])*vscale;
		tgp.vel[j] += hubble*tgp.pos[j];
		}
	    tgp.mass *= mscale;
	    tgp.rho *= mscale/(rscale*rscale*rscale);
	    tgp.hsmooth *= rscale;
	    tgp.phi *= vscale*vscale;
	    write_tipsy_xdr_gas(&xdrsout,&tgp);
	    }
	for(i = 0; i < th.ndark; i++) {
	    read_tipsy_xdr_dark(&xdrsin,&tdp);
	    for(j = 0; j < 3; j++) {
		tdp.pos[j] = (tdp.pos[j]-rcen[j])*rscale;
		tdp.vel[j] = (tdp.vel[j]-vcen[j])*vscale;
		tdp.vel[j] += hubble*tdp.pos[j];
		}
	    tdp.mass *= mscale;
	    tdp.eps *= rscale;
	    tdp.phi *= vscale*vscale;
	    write_tipsy_xdr_dark(&xdrsout,&tdp);
	    }
	for(i = 0; i < th.nstar; i++) {
	    read_tipsy_xdr_star(&xdrsin,&tsp);
	    for(j = 0; j < 3; j++) {
		tsp.pos[j] = (tsp.pos[j]-rcen[j])*rscale;
		tsp.vel[j] = (tsp.vel[j]-vcen[j])*vscale;
		tsp.vel[j] += hubble*tsp.pos[j];
		}
	    tsp.mass *= mscale;
	    tsp.eps *= rscale;
	    tsp.phi *= vscale*vscale;
	    write_tipsy_xdr_star(&xdrsout,&tsp);
	    }
	}
    else if (positionprecision == 1) {
	for(i = 0; i < th.ngas; i++) {
	    read_tipsy_xdr_gas_dpp(&xdrsin,&tgpdpp);
	    for(j = 0; j < 3; j++) {
		tgpdpp.pos[j] = (tgpdpp.pos[j]-rcen[j])*rscale;
		tgpdpp.vel[j] = (tgpdpp.vel[j]-vcen[j])*vscale;
		tgpdpp.vel[j] += hubble*tgpdpp.pos[j];
		}
	    tgpdpp.mass *= mscale;
	    tgpdpp.rho *= mscale/(rscale*rscale*rscale);
	    tgpdpp.hsmooth *= rscale;
	    tgpdpp.phi *= vscale*vscale;
	    write_tipsy_xdr_gas_dpp(&xdrsout,&tgpdpp);
	    }
	for(i = 0; i < th.ndark; i++) {
	    read_tipsy_xdr_dark_dpp(&xdrsin,&tdpdpp);
	    for(j = 0; j < 3; j++) {
		tdpdpp.pos[j] = (tdpdpp.pos[j]-rcen[j])*rscale;
		tdpdpp.vel[j] = (tdpdpp.vel[j]-vcen[j])*vscale;
		tdpdpp.vel[j] += hubble*tdpdpp.pos[j];
		}
	    tdpdpp.mass *= mscale;
	    tdpdpp.eps *= rscale;
	    tdpdpp.phi *= vscale*vscale;
	    write_tipsy_xdr_dark_dpp(&xdrsout,&tdpdpp);
	    }
	for(i = 0; i < th.nstar; i++) {
	    read_tipsy_xdr_star_dpp(&xdrsin,&tspdpp);
	    for(j = 0; j < 3; j++) {
		tspdpp.pos[j] = (tspdpp.pos[j]-rcen[j])*rscale;
		tspdpp.vel[j] = (tspdpp.vel[j]-vcen[j])*vscale;
		tspdpp.vel[j] += hubble*tspdpp.pos[j];
		}
	    tspdpp.mass *= mscale;
	    tspdpp.eps *= rscale;
	    tspdpp.phi *= vscale*vscale;
	    write_tipsy_xdr_star_dpp(&xdrsout,&tspdpp);
	    }
	}
    /*
    ** Clean up and write some output
    */
    xdr_destroy(&xdrsin);
    xdr_destroy(&xdrsout);
    if (verboselevel >= 1) {
	fprintf(stderr,"Used values:\n");
	fprintf(stderr,"rcen   : (%+.6e,%+.6e,%+.6e) LUold\n",rcen[0],rcen[1],rcen[2]);
	fprintf(stderr,"vcen   : (%+.6e,%+.6e,%+.6e) VUold\n",vcen[0],vcen[1],vcen[2]);
	fprintf(stderr,"mscale : %.6e MUnew/MUold\n",mscale);
	fprintf(stderr,"rscale : %.6e LUnew/LUold\n",rscale);
	fprintf(stderr,"vscale : %.6e VUnew/VUold\n",vscale);
	fprintf(stderr,"Hubble : %.6e VUnew/LUnew\n",hubble);
	}
    if (verboselevel >= 0) {
	fprintf(stderr,"Time: %g Ntotal: %u Ngas: %u Ndark: %u Nstar: %u\n",
		th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
	}
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program makes a coordinate transformation from comoving to physical coordinates.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"r_new = (r_old-r_cen)*rscale (position)\n");
    fprintf(stderr,"v_new = (v_old-v_cen)*vscale + Hubble*r_new (velocity)\n");
    fprintf(stderr,"m_new = m_old*mscale (mass)\n");
    fprintf(stderr,"l_new = l_old*rscale (length)\n");
    fprintf(stderr,"E_new = E_old*vscale^2 (energy/mass)\n");
    fprintf(stderr,"Temperature, metals and formation time are left unchanged.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp            : set this flag if input and output files have single precision positions (default)\n");
    fprintf(stderr,"-dpp            : set this flag if input and output files have double precision positions\n");
    fprintf(stderr,"-rxcen <value>  : x-coordinate of centre [LUold] (default: 0 LUold)\n");
    fprintf(stderr,"-rycen <value>  : y-coordinate of centre [LUold] (default: 0 LUold)\n");
    fprintf(stderr,"-rzcen <value>  : z-coordinate of centre [LUold] (default: 0 LUold)\n");
    fprintf(stderr,"-vxcen <value>  : x-velocity of centre [VUold] (default: 0 VUold)\n");
    fprintf(stderr,"-vycen <value>  : y-velocity of centre [VUold] (default: 0 VUold)\n");
    fprintf(stderr,"-vzcen <value>  : z-velocity of centre [VUold] (default: 0 VUold)\n");
    fprintf(stderr,"-mscale <value> : mass scale factor [MUnew/MUold] (default: 1 MUnew/MUold)\n");
    fprintf(stderr,"-rscale <value> : coordinate scale factor [LUnew/LUold] (default: 1 LUnew/LUold)\n");
    fprintf(stderr,"-vscale <value> : velocity scale factor [VUnew/VUold] (default: 1 VUnew/VUold)\n");
    fprintf(stderr,"-Hubble <value> : Hubble parameter [1/TUnew] (default: 0 1/TUnew; special value: sqrt_8pi_3)\n");
    fprintf(stderr,"-v              : more informative output to screen\n");
    fprintf(stderr,"< <name>        : input file in tipsy XDR format\n");
    fprintf(stderr,"> <name>        : output file in tipsy XDR format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

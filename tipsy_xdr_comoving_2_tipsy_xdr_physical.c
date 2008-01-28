/* 
** tscom2tsphy.c
**
** written by Marcel Zemp
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
    int positionprecision;
    double rcen[3], vcen[3];
    double mscale, rscale, vscale, hubble;
    TIPSY_HEADER th;
    GAS_PARTICLE gp;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
    GAS_PARTICLE_DPP gpdpp;
    DARK_PARTICLE_DPP dpdpp;
    STAR_PARTICLE_DPP spdpp;
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
    xdrstdio_create(&xdrsin,stdin,XDR_DECODE);
    xdrstdio_create(&xdrsout,stdout,XDR_ENCODE);
    read_tipsy_standard_header(&xdrsin,&th);
    write_tipsy_standard_header(&xdrsout,&th);
    /*
    ** Add shifts to particles and write them out
    */
    if (positionprecision == 0) {
	for(i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas(&xdrsin,&gp);
	    for(j = 0; j < 3; j++) {
		gp.pos[j] = (gp.pos[j]-rcen[j])*rscale;
		gp.vel[j] = (gp.vel[j]-vcen[j])*vscale;
		gp.vel[j] += hubble*gp.pos[j];
		}
	    gp.mass *= mscale;
	    gp.rho *= mscale/(rscale*rscale*rscale);
	    gp.hsmooth *= rscale;
	    gp.phi *= vscale*vscale;
	    write_tipsy_standard_gas(&xdrsout,&gp);
	    }
	for(i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark(&xdrsin,&dp);
	    for(j = 0; j < 3; j++) {
		dp.pos[j] = (dp.pos[j]-rcen[j])*rscale;
		dp.vel[j] = (dp.vel[j]-vcen[j])*vscale;
		dp.vel[j] += hubble*dp.pos[j];
		}
	    dp.mass *= mscale;
	    dp.eps *= rscale;
	    dp.phi *= vscale*vscale;
	    write_tipsy_standard_dark(&xdrsout,&dp);
	    }
	for(i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star(&xdrsin,&sp);
	    for(j = 0; j < 3; j++) {
		sp.pos[j] = (sp.pos[j]-rcen[j])*rscale;
		sp.vel[j] = (sp.vel[j]-vcen[j])*vscale;
		sp.vel[j] += hubble*sp.pos[j];
		}
	    sp.mass *= mscale;
	    sp.eps *= rscale;
	    sp.phi *= vscale*vscale;
	    write_tipsy_standard_star(&xdrsout,&sp);
	    }
	}
    else if (positionprecision == 1) {
	for(i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas_dpp(&xdrsin,&gpdpp);
	    for(j = 0; j < 3; j++) {
		gpdpp.pos[j] = (gpdpp.pos[j]-rcen[j])*rscale;
		gpdpp.vel[j] = (gpdpp.vel[j]-vcen[j])*vscale;
		gpdpp.vel[j] += hubble*gpdpp.pos[j];
		}
	    gpdpp.mass *= mscale;
	    gpdpp.rho *= mscale/(rscale*rscale*rscale);
	    gpdpp.hsmooth *= rscale;
	    gpdpp.phi *= vscale*vscale;
	    write_tipsy_standard_gas_dpp(&xdrsout,&gpdpp);
	    }
	for(i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark_dpp(&xdrsin,&dpdpp);
	    for(j = 0; j < 3; j++) {
		dpdpp.pos[j] = (dpdpp.pos[j]-rcen[j])*rscale;
		dpdpp.vel[j] = (dpdpp.vel[j]-vcen[j])*vscale;
		dpdpp.vel[j] += hubble*dpdpp.pos[j];
		}
	    dpdpp.mass *= mscale;
	    dpdpp.eps *= rscale;
	    dpdpp.phi *= vscale*vscale;
	    write_tipsy_standard_dark_dpp(&xdrsout,&dpdpp);
	    }
	for(i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star_dpp(&xdrsin,&spdpp);
	    for(j = 0; j < 3; j++) {
		spdpp.pos[j] = (spdpp.pos[j]-rcen[j])*rscale;
		spdpp.vel[j] = (spdpp.vel[j]-vcen[j])*vscale;
		spdpp.vel[j] += hubble*spdpp.pos[j];
		}
	    spdpp.mass *= mscale;
	    spdpp.eps *= rscale;
	    spdpp.phi *= vscale*vscale;
	    write_tipsy_standard_star_dpp(&xdrsout,&spdpp);
	    }
	}
    /*
    ** Clean up and write some output
    */
    xdr_destroy(&xdrsin);
    xdr_destroy(&xdrsout);
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
    fprintf(stderr,"Used values:\n");
    fprintf(stderr,"rcen   = (%+.6e,%+.6e,%+.6e) LUold\n",rcen[0],rcen[1],rcen[2]);
    fprintf(stderr,"vcen   = (%+.6e,%+.6e,%+.6e) VUold\n",vcen[0],vcen[1],vcen[2]);
    fprintf(stderr,"mscale = %.6e MUnew/MUold\n",mscale);
    fprintf(stderr,"rscale = %.6e LUnew/LUold\n",rscale);
    fprintf(stderr,"vscale = %.6e VUnew/VUold\n",vscale);
    fprintf(stderr,"Hubble = %.6e VUnew/LUnew\n",hubble);
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
    fprintf(stderr,"-spp            : set this flag if input and output file have single precision positions (default)\n");
    fprintf(stderr,"-dpp            : set this flag if input and output file have double precision positions\n");
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
    fprintf(stderr,"< <name>        : input file in tipsy standard binary format\n");
    fprintf(stderr,"> <name>        : output file in tipsy standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

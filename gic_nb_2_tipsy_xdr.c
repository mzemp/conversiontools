/* 
** gicb2ts.c
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <iof.h>
#include <gic_reader.h>

void usage(void);

int main(int argc, char **argv) {

    int i = 0, j = 0, k = 0;
    int positionprecision, verboselevel, multires, particletype;
    long Ntot, Nlev, index;
    int Nrec, Npage, N1D, Nlow, Lmax;
    GIC_REAL *PosXRec, *PosYRec, *PosZRec;
    GIC_REAL *VelXRec, *VelYRec, *VelZRec;
    double posscalefac, velscalefac, particlemass, particlesoftening;
    double dr[3], dx, aBegin, h100, H_Tipsy, TU_Tipsy, LU_Tipsy, VU_Tipsy, MU_Tipsy, VelConvertFac, rhocrit;
    char PosFileName[256], VelFileName[256];
    TIPSY_HEADER th;
    GAS_PARTICLE gp;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
    GAS_PARTICLE_DPP gpdpp;
    DARK_PARTICLE_DPP dpdpp;
    STAR_PARTICLE_DPP spdpp;
    struct gicFile PosXFile, PosYFile, PosZFile;
    struct gicFile VelXFile, VelYFile, VelZFile;
    struct gicManifest PosManifest, VelManifest;
    struct gicFileHeader PosHeader, VelHeader;
    struct gicLevelHeader *PosLevelHeader = NULL, *VelLevelHeader = NULL;
    XDR xdrs;

    /*
    ** Set some default values
    */

    positionprecision = 0;
    verboselevel = 0;
    multires = 0;
    particletype = 1;
    dr[0] = -0.5;
    dr[1] = -0.5;
    dr[2] = -0.5;
    particlesoftening = -1;
    posscalefac = -1;
    velscalefac = -1;
    H_Tipsy = sqrt(8*M_PI/3); /* TU_Tipsy^-1 */
    VelConvertFac = 1.0227121651152353693;
    rhocrit = 277.53662719; /* h^2 Mo kpc^-3 */

    /*
    ** Read in options
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
	else if (strcmp(argv[i],"-sr") == 0) {
            multires = 0;
            i++;
            }
	else if (strcmp(argv[i],"-mr") == 0) {
            multires = 1;
            i++;
            }
	else if (strcmp(argv[i],"-gas") == 0) {
            particletype = 0;
            i++;
            }
	else if (strcmp(argv[i],"-dark") == 0) {
            particletype = 1;
            i++;
            }
	else if (strcmp(argv[i],"-star") == 0) {
            particletype = 2;
            i++;
            }
        else if (strcmp(argv[i],"-drx") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            dr[0] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-dry") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            dr[1] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-drz") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            dr[2] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-soft") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            particlesoftening = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-mass") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            particlemass = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-posfac") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            posscalefac = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-velfac") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            velscalefac = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-pos") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            strcpy(PosFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-vel") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            strcpy(VelFileName,argv[i]);
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

    /*
    ** Read in manifests & headers of the two files. 
    ** Attention you need to read through since Nrec and WrongOrder are set by this process.
    */

    PosXFile.File = fopen(PosFileName,"r");
    assert(PosXFile.File != NULL);
    PosYFile.File = fopen(PosFileName,"r");
    assert(PosYFile.File != NULL);
    PosZFile.File = fopen(PosFileName,"r");
    assert(PosZFile.File != NULL);

    VelXFile.File = fopen(VelFileName,"r");
    assert(VelXFile.File != NULL);
    VelYFile.File = fopen(VelFileName,"r");
    assert(VelYFile.File != NULL);
    VelZFile.File = fopen(VelFileName,"r");
    assert(VelZFile.File != NULL);

    assert(gicReadManifest(&PosXFile,&PosManifest) == 0);
    assert(gicReadManifest(&PosYFile,&PosManifest) == 0);
    assert(gicReadManifest(&PosZFile,&PosManifest) == 0);

    assert(gicReadManifest(&VelXFile,&VelManifest) == 0);
    assert(gicReadManifest(&VelYFile,&VelManifest) == 0);
    assert(gicReadManifest(&VelZFile,&VelManifest) == 0);

    assert(gicReadFileHeader(&PosXFile,&PosHeader) == 0);
    assert(gicReadFileHeader(&PosYFile,&PosHeader) == 0);
    assert(gicReadFileHeader(&PosZFile,&PosHeader) == 0);

    assert(gicReadFileHeader(&VelXFile,&VelHeader) == 0);
    assert(gicReadFileHeader(&VelYFile,&VelHeader) == 0);
    assert(gicReadFileHeader(&VelZFile,&VelHeader) == 0);

    /*
    ** Do some consitency tests, set some variables and tipsy header
    */

    assert(PosXFile.Nrec == VelXFile.Nrec);
    assert(PosManifest.OmegaB == VelManifest.OmegaB);
    assert(PosManifest.OmegaX == VelManifest.OmegaX);
    assert(PosManifest.dx == VelManifest.dx);
    assert(PosManifest.h100 == VelManifest.h100);
    assert(PosHeader.aBegin == VelHeader.aBegin);
    assert(PosHeader.Ntot == VelHeader.Ntot);
    assert(PosHeader.Lmax == VelHeader.Lmax);
    assert(PosHeader.dims[0] == VelHeader.dims[0]);
    assert(PosHeader.dims[1] == VelHeader.dims[1]);
    assert(PosHeader.dims[2] == VelHeader.dims[2]);

    Ntot = PosHeader.Ntot;
    Nlow = PosHeader.dims[0]*PosHeader.dims[1]*PosHeader.dims[2];
    Nrec = PosXFile.Nrec;
    N1D = PosHeader.dims[0];
    dx = PosManifest.dx;
    aBegin = PosHeader.aBegin;
    h100 = PosManifest.h100;
    Lmax = PosHeader.Lmax;

    th.time = aBegin;
    th.ntotal = Ntot;
    th.ndim = 3;
    if (particletype == 0) {
	th.ngas = Ntot; 
	th.ndark = 0;
	th.nstar = 0;
	}
    else if (particletype == 1) {
	th.ngas = 0;
	th.ndark = Ntot;
	th.nstar = 0;
	}
    else {
	th.ngas = 0;
	th.ndark = 0;
	th.nstar = Ntot;
	}

    /*
    ** Allocate level headers and read-in records
    */

    PosLevelHeader = malloc((Lmax+1)*sizeof(struct gicLevelHeader));
    VelLevelHeader = malloc((Lmax+1)*sizeof(struct gicLevelHeader));

    PosXRec = malloc(Nrec*sizeof(GIC_REAL));
    PosYRec = malloc(Nrec*sizeof(GIC_REAL));
    PosZRec = malloc(Nrec*sizeof(GIC_REAL));

    VelXRec = malloc(Nrec*sizeof(GIC_REAL));
    VelYRec = malloc(Nrec*sizeof(GIC_REAL));
    VelZRec = malloc(Nrec*sizeof(GIC_REAL));

    /*
    ** Calculate units and scaling factors
    */

    TU_Tipsy = (H_Tipsy*1000)/(h100*100); /* kpc*km^-1*s */
    LU_Tipsy = N1D*dx*1000/h100; /* kpc */
    VU_Tipsy = LU_Tipsy/TU_Tipsy; /* km s^-1 */
    MU_Tipsy = rhocrit*h100*h100*LU_Tipsy*LU_Tipsy*LU_Tipsy; /* Mo */

    if (posscalefac < 0) posscalefac = 1/(N1D*dx);
    if (velscalefac < 0) velscalefac = 1/(aBegin*VU_Tipsy);

    /*
    ** Get output file ready
    */

    xdrstdio_create(&xdrs,stdout,XDR_ENCODE);
    write_tipsy_standard_header(&xdrs,&th);

    /*
    ** Read in data
    */

    if (multires == 0 ) { /* single resolution */

	/*
	** Set pointers to streams to the correct positions
	*/

	Npage = (Ntot+Nrec-1)/Nrec;
	for (i = 0; i < Npage; i++) {
	    assert(gicSkipFortranRecordReal(&PosYFile) == 0);
	    assert(gicSkipFortranRecordReal(&PosZFile) == 0);
	    assert(gicSkipFortranRecordReal(&VelYFile) == 0);
	    assert(gicSkipFortranRecordReal(&VelZFile) == 0);
	    }
	for (i = 0; i < Npage; i++) {
	    assert(gicSkipFortranRecordReal(&PosZFile) == 0);
	    assert(gicSkipFortranRecordReal(&VelZFile) == 0);
	    }

	/*
	** Set particle mass and softening
	*/

	if (particlesoftening < 0) particlesoftening = 1.0/(N1D*20.0);
	if (particlemass < 0) particlemass = (PosManifest.OmegaB+PosManifest.OmegaX)/Ntot;

	/*
	** Process the data
	*/

	index = 0;
	for (i = 0; i < Npage; i++) {
	    assert(gicReadFortranRecordReal(&PosXFile,PosXRec) == 0);
	    assert(gicReadFortranRecordReal(&PosYFile,PosYRec) == 0);
	    assert(gicReadFortranRecordReal(&PosZFile,PosZRec) == 0);
	    assert(gicReadFortranRecordReal(&VelXFile,VelXRec) == 0);
	    assert(gicReadFortranRecordReal(&VelYFile,VelYRec) == 0);
	    assert(gicReadFortranRecordReal(&VelZFile,VelZRec) == 0);
	    for (j = 0; j < Nrec; j++) {
		if (positionprecision == 0) {
		    if (particletype == 0) {
			gp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			gp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			gp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			gp.vel[0] = VelXRec[j]*velscalefac;
			gp.vel[1] = VelYRec[j]*velscalefac;
			gp.vel[2] = VelZRec[j]*velscalefac;
			gp.mass = particlemass;
			gp.rho = 0;
			gp.temp = 0;
			gp.hsmooth = particlesoftening;
			gp.metals = 0;
			gp.phi = 0;
			if (index < Ntot) write_tipsy_standard_gas(&xdrs,&gp);
			}
		    else if (particletype == 1) {
			dp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			dp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			dp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			dp.vel[0] = VelXRec[j]*velscalefac;
			dp.vel[1] = VelYRec[j]*velscalefac;
			dp.vel[2] = VelZRec[j]*velscalefac;
			dp.mass = particlemass;
			dp.eps = particlesoftening;
			dp.phi = 0;
			if (index < Ntot) write_tipsy_standard_dark(&xdrs,&dp);
			}
		    else {
			sp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			sp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			sp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			sp.vel[0] = VelXRec[j]*velscalefac;
			sp.vel[1] = VelYRec[j]*velscalefac;
			sp.vel[2] = VelZRec[j]*velscalefac;
			sp.mass = particlemass;
			sp.metals = 0;
			sp.tform = 0;
			sp.eps = particlesoftening;
			sp.phi = 0;
			if (index < Ntot) write_tipsy_standard_star(&xdrs,&sp);
			}
		    }
		else {
		    if (particletype == 0) {
			gpdpp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			gpdpp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			gpdpp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			gpdpp.vel[0] = VelXRec[j]*velscalefac;
			gpdpp.vel[1] = VelYRec[j]*velscalefac;
			gpdpp.vel[2] = VelZRec[j]*velscalefac;
			gpdpp.mass = particlemass;
			gpdpp.rho = 0;
			gpdpp.temp = 0;
			gpdpp.hsmooth = particlesoftening;
			gpdpp.metals = 0;
			gpdpp.phi = 0;
			if (index < Ntot) write_tipsy_standard_gas_dpp(&xdrs,&gpdpp);
			}
		    else if (particletype == 1) {
			dpdpp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			dpdpp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			dpdpp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			dpdpp.vel[0] = VelXRec[j]*velscalefac;
			dpdpp.vel[1] = VelYRec[j]*velscalefac;
			dpdpp.vel[2] = VelZRec[j]*velscalefac;
			dpdpp.mass = particlemass;
			dpdpp.eps = particlesoftening;
			dpdpp.phi = 0;
			if (index < Ntot) write_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
			}
		    else {
			spdpp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			spdpp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			spdpp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			spdpp.vel[0] = VelXRec[j]*velscalefac;
			spdpp.vel[1] = VelYRec[j]*velscalefac;
			spdpp.vel[2] = VelZRec[j]*velscalefac;
			spdpp.mass = particlemass;
			spdpp.metals = 0;
			spdpp.tform = 0;
			spdpp.eps = particlesoftening;
			spdpp.phi = 0;
			if (index < Ntot) write_tipsy_standard_star_dpp(&xdrs,&spdpp);
			}
		    }
		index++;
		}
	    }
	}
    else { /* multi-resolution */

	/*
	** Some info
	*/

	if (verboselevel >= 1) {
	    fprintf(stderr,"There are %d refinement levels:\n\n",Lmax+1);
	    }

	/*
	** Skip refinement mask array
	*/

	Npage = (Nlow+Nrec-1)/Nrec;
	for (i = 0; i < Npage; i++) {
	    assert(gicSkipFortranRecordInteger(&PosXFile) == 0);
	    assert(gicSkipFortranRecordInteger(&PosYFile) == 0);
	    assert(gicSkipFortranRecordInteger(&PosZFile) == 0);
	    assert(gicSkipFortranRecordInteger(&VelXFile) == 0);
	    assert(gicSkipFortranRecordInteger(&VelYFile) == 0);
	    assert(gicSkipFortranRecordInteger(&VelZFile) == 0);
	    }

	/*
	** Read levels
	*/

	for (k = 0; k <= Lmax; k++) {

	    /*
	    ** Read level headers
	    */

	    assert(gicReadLevelHeader(&PosXFile,&PosLevelHeader[k]) == 0);
	    assert(gicReadLevelHeader(&PosYFile,&PosLevelHeader[k]) == 0);
	    assert(gicReadLevelHeader(&PosZFile,&PosLevelHeader[k]) == 0);
	    assert(gicReadLevelHeader(&VelXFile,&VelLevelHeader[k]) == 0);
	    assert(gicReadLevelHeader(&VelYFile,&VelLevelHeader[k]) == 0);
	    assert(gicReadLevelHeader(&VelZFile,&VelLevelHeader[k]) == 0);

	    /*
	    ** Set pointers to streams to the correct positions
	    */
	    
	    assert(PosLevelHeader[k].Nlev == VelLevelHeader[k].Nlev);
	    Nlev = PosLevelHeader[k].Nlev;
	    Npage = (Nlev+Nrec-1)/Nrec;
	    for (i = 0; i < Npage; i++) {
		assert(gicSkipFortranRecordReal(&PosYFile) == 0);
		assert(gicSkipFortranRecordReal(&PosZFile) == 0);
		assert(gicSkipFortranRecordReal(&VelYFile) == 0);
		assert(gicSkipFortranRecordReal(&VelZFile) == 0);
		}
	    for (i = 0; i < Npage; i++) {
		assert(gicSkipFortranRecordReal(&PosZFile) == 0);
		assert(gicSkipFortranRecordReal(&VelZFile) == 0);
		}

	    /*
	    ** Set particle softening
	    */

	    if (k == 0) { /* top level softening */
		if (particlesoftening < 0) particlesoftening = 1.0/(N1D*20.0);
		}
	    else {
		particlesoftening /= 2.0;
		}
	    assert(PosLevelHeader[k].Mlev == VelLevelHeader[k].Mlev);
	    particlemass = PosLevelHeader[k].Mlev/MU_Tipsy;

	    /*
	    ** Process the data
	    */
	    
	    index = 0;
	    for (i = 0; i < Npage; i++) {
		assert(gicReadFortranRecordReal(&PosXFile,PosXRec) == 0);
		assert(gicReadFortranRecordReal(&PosYFile,PosYRec) == 0);
		assert(gicReadFortranRecordReal(&PosZFile,PosZRec) == 0);
		assert(gicReadFortranRecordReal(&VelXFile,VelXRec) == 0);
		assert(gicReadFortranRecordReal(&VelYFile,VelYRec) == 0);
		assert(gicReadFortranRecordReal(&VelZFile,VelZRec) == 0);
		for (j = 0; j < Nrec; j++) {
		    if (positionprecision == 0) {
			if (particletype == 0) {
			    gp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			    gp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			    gp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			    gp.vel[0] = VelXRec[j]*velscalefac;
			    gp.vel[1] = VelYRec[j]*velscalefac;
			    gp.vel[2] = VelZRec[j]*velscalefac;
			    gp.mass = particlemass;
			    gp.rho = 0;
			    gp.temp = 0;
			    gp.hsmooth = particlesoftening;
			    gp.metals = 0;
			    gp.phi = 0;
			    if (index < Nlev) write_tipsy_standard_gas(&xdrs,&gp);
			    }
			else if (particletype == 1) {
			    dp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			    dp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			    dp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			    dp.vel[0] = VelXRec[j]*velscalefac;
			    dp.vel[1] = VelYRec[j]*velscalefac;
			    dp.vel[2] = VelZRec[j]*velscalefac;
			    dp.mass = particlemass;
			    dp.eps = particlesoftening;
			    dp.phi = 0;
			    if (index < Nlev) write_tipsy_standard_dark(&xdrs,&dp);
			    }
			else {
			    sp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			    sp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			    sp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			    sp.vel[0] = VelXRec[j]*velscalefac;
			    sp.vel[1] = VelYRec[j]*velscalefac;
			    sp.vel[2] = VelZRec[j]*velscalefac;
			    sp.mass = particlemass;
			    sp.metals = 0;
			    sp.tform = 0;
			    sp.eps = particlesoftening;
			    sp.phi = 0;
			    if (index < Nlev) write_tipsy_standard_star(&xdrs,&sp);
			    }
			}
		    else {
			if (particletype == 0) {
			    gpdpp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			    gpdpp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			    gpdpp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			    gpdpp.vel[0] = VelXRec[j]*velscalefac;
			    gpdpp.vel[1] = VelYRec[j]*velscalefac;
			    gpdpp.vel[2] = VelZRec[j]*velscalefac;
			    gpdpp.mass = particlemass;
			    gpdpp.rho = 0;
			    gpdpp.temp = 0;
			    gpdpp.hsmooth = particlesoftening;
			    gpdpp.metals = 0;
			    gpdpp.phi = 0;
			    if (index < Nlev) write_tipsy_standard_gas_dpp(&xdrs,&gpdpp);
			    }
			else if (particletype == 1) {
			    dpdpp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			    dpdpp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			    dpdpp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			    dpdpp.vel[0] = VelXRec[j]*velscalefac;
			    dpdpp.vel[1] = VelYRec[j]*velscalefac;
			    dpdpp.vel[2] = VelZRec[j]*velscalefac;
			    dpdpp.mass = particlemass;
			    dpdpp.eps = particlesoftening;
			    dpdpp.phi = 0;
			    if (index < Nlev) write_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
			    }
			else {
			    spdpp.pos[0] = PosXRec[j]*posscalefac + dr[0];
			    spdpp.pos[1] = PosYRec[j]*posscalefac + dr[1];
			    spdpp.pos[2] = PosZRec[j]*posscalefac + dr[2];
			    spdpp.vel[0] = VelXRec[j]*velscalefac;
			    spdpp.vel[1] = VelYRec[j]*velscalefac;
			    spdpp.vel[2] = VelZRec[j]*velscalefac;
			    spdpp.mass = particlemass;
			    spdpp.metals = 0;
			    spdpp.tform = 0;
			    spdpp.eps = particlesoftening;
			    spdpp.phi = 0;
			    if (index < Nlev) write_tipsy_standard_star_dpp(&xdrs,&spdpp);
			    }
			}
		    index++;
		    }
		}
	    
	    /*
	    ** Get streams ready for next level
	    */
	    
	    for (i = 0; i < Npage; i++) {
		assert(gicSkipFortranRecordReal(&PosXFile) == 0);
		assert(gicSkipFortranRecordReal(&PosYFile) == 0);
		assert(gicSkipFortranRecordReal(&VelXFile) == 0);
		assert(gicSkipFortranRecordReal(&VelYFile) == 0);
		}
	    for (i = 0; i < Npage; i++) {
		assert(gicSkipFortranRecordReal(&PosXFile) == 0);
		assert(gicSkipFortranRecordReal(&VelXFile) == 0);
		}

	    /*
	    ** Some output about levels
	    */

	    if (verboselevel >= 1) {
		fprintf(stderr,"L %d Lmax %d Nlev %ld soft %g mass %g\n",PosLevelHeader[k].L,PosLevelHeader[k].Lmax,Nlev,particlesoftening,particlemass);
		if (k == Lmax) fprintf(stderr,"\n");
		}
	    }
	}

    xdr_destroy(&xdrs);

    /*
    ** Write out some additional stuff depending on verbose level
    */

    if (verboselevel >= 1) {
	fprintf(stderr,"Parameters from general initial conditions file:\n\n");
	fprintf(stderr,"Name     : %s\n",PosManifest.name);
	fprintf(stderr,"OmegaB   : %.6e\n",PosManifest.OmegaB);
	fprintf(stderr,"OmegaDM  : %.6e\n",PosManifest.OmegaX);
	fprintf(stderr,"OmegaL   : %.6e\n",PosManifest.OmegaL);
	fprintf(stderr,"OmegaN   : %.6e\n",PosManifest.OmegaN);
	fprintf(stderr,"h100     : %.6e\n",h100);
	fprintf(stderr,"Deltax   : %.6e chimp\n",dx);
	fprintf(stderr,"ns       : %.6e\n",PosManifest.ns);
	fprintf(stderr,"sigma8   : %.6e\n",PosManifest.s8);
	fprintf(stderr,"kp       : %.6e\n",PosManifest.kp);
	fprintf(stderr,"aBegin   : %.6e\n",aBegin);
	fprintf(stderr,"DeltaDC  : %.6e\n",PosHeader.DeltaDC);
	fprintf(stderr,"NX       : %d\n",PosHeader.dims[0]);
	fprintf(stderr,"NY       : %d\n",PosHeader.dims[1]);
	fprintf(stderr,"NZ       : %d\n",PosHeader.dims[2]);
	fprintf(stderr,"Seed     : %d\n",PosHeader.seed);
	fprintf(stderr,"Nrec     : %d\n",Nrec);
	fprintf(stderr,"Ntot     : %ld\n",Ntot);
	fprintf(stderr,"Lmax     : %d\n\n",Lmax);
	fprintf(stderr,"Used values:\n\n");
	fprintf(stderr,"drx  : %.6e LU\n",dr[0]);
	fprintf(stderr,"dry  : %.6e LU\n",dr[1]);
	fprintf(stderr,"drz  : %.6e LU\n",dr[2]);
	fprintf(stderr,"soft : %.6e LU\n\n",particlesoftening*pow(2,k));
	fprintf(stderr,"Resulting internal tipsy units:\n\n");
	fprintf(stderr,"LU : %.6e kpc\n",LU_Tipsy);
	fprintf(stderr,"TU : %.6e Gyr\n",TU_Tipsy/VelConvertFac);
	fprintf(stderr,"VU : %.6e km s^-1 = %.6e kpc Gyr^-1\n",VU_Tipsy,VU_Tipsy*VelConvertFac);
	fprintf(stderr,"MU : %.6e Mo\n\n",MU_Tipsy);
	}
    if (verboselevel >= 0) {
	fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	       th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
	}

    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts general initial conditions format to tipsy standard binary format.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp            : set this flag if input and output files have single precision positions (default)\n");
    fprintf(stderr,"-dpp            : set this flag if input and output files have double precision positions\n");
    fprintf(stderr,"-sr             : set this flag if input only has single resolution particles (default)\n");
    fprintf(stderr,"-mr             : set this flag if input has multi resolution particles\n");
    fprintf(stderr,"-gas            : set this flag if input particles are gas particles\n");
    fprintf(stderr,"-dark           : set this flag if input particles are dark particles (default)\n");
    fprintf(stderr,"-star           : set this flag if input particles are star particles\n");
    fprintf(stderr,"-drx <value>    : shift along x-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-dry <value>    : shift along y-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-drz <value>    : shift along z-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-soft <value>   : softening length of particles [LU] (default: 1/20 mean particle separation => 1/(N^{-3}*20) LU)\n");
    fprintf(stderr,"-mass <value>   : mass of particles [MU] (default: OmegaM/Ntot [sr] or from file [mr])\n");
    fprintf(stderr,"-posfac <value> : position scale factor (default: 1/[N1D*dx] where dx is the cell size)\n");
    fprintf(stderr,"-velfac <value> : velocity scale factor (default: 1/[a*VU_Tipsy] where a is the scale factor)\n");
    fprintf(stderr,"-v              : more informative output to screen\n");
    fprintf(stderr,"-pos <name>     : positions input file gic binary format\n");
    fprintf(stderr,"-vel <name>     : velocities input file gic binary format\n");
    fprintf(stderr,"> <name>        : output file in tipsy standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

/* 
** gic_nb_2_tipsy_xdr.c
**
** Written by Marcel Zemp
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
#include <gicreader.h>

void usage(void);

int main(int argc, char **argv) {

	int i = 0, j = 0, k = 0;
	int positionprecision, verboselevel, multires, particletype, mrmassfromfile;
	long Ntot, Nlev, index;
	int Nrec, Npage, N1Dlow, Nlow, Lmax;
	GIC_REAL *PosXRec, *PosYRec, *PosZRec;
	GIC_REAL *VelXRec, *VelYRec, *VelZRec;
	double LUsf, VUsf, MUsf;
	double b2dmscalefac, csvelscalefac, refinementstep;
	double toplevelmass, toplevelsoftening;
	double mass, soft;
	double dr[3], dx, acurrent, VelConvertFac, rhocrit, LBox, softfac;
	double OmegaM0, OmegaDM0, OmegaB0, OmegaL0, OmegaK0, OmegaR0, h100;
	double H_TIPSY_DEFAULT, TU_TIPSY_DEFAULT, LU_TIPSY_DEFAULT, VU_TIPSY_DEFAULT, MU_TIPSY_DEFAULT;
	double TU_TIPSY, LU_TIPSY, VU_TIPSY, MU_TIPSY;
	double TU_GIC, LU_GIC, VU_GIC, MU_GIC;
	char PosFileName[256], VelFileName[256];
	TIPSY_HEADER th;
	TIPSY_GAS_PARTICLE tgp;
	TIPSY_DARK_PARTICLE tdp;
	TIPSY_STAR_PARTICLE tsp;
	TIPSY_GAS_PARTICLE_DPP tgpdpp;
	TIPSY_DARK_PARTICLE_DPP tdpdpp;
	TIPSY_STAR_PARTICLE_DPP tspdpp;
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
	mrmassfromfile = 0;
	dr[0] = -0.5;
	dr[1] = -0.5;
	dr[2] = -0.5;
	refinementstep = 2;
	toplevelsoftening = -2;
	toplevelmass = -2;
	LUsf = -2;
	VUsf = -2;
	MUsf = -2;
	softfac = 50;
	b2dmscalefac = -1;
	H_TIPSY_DEFAULT = sqrt(8*M_PI/3); /* TU_TIPSY^-1 */
	VelConvertFac = ConversionFactors.km_per_s_2_kpc_per_Gyr;
	rhocrit = PhysicalConstants.rho_crit_cosmology; /* h^2 Mo kpc^-3 */

	/*
	** Read in options
	*/

	i = 1;
	while (i < argc) {
		if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
			usage();
			}
		if (strcmp(argv[i],"-version") == 0) {
			fprintf(stderr,"%s (%s)\n",NAME,VERSION);
			exit(1);
			}
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
			if (i >= argc) usage();
			dr[0] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-dry") == 0) {
			i++;
			if (i >= argc) usage();
			dr[1] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-drz") == 0) {
			i++;
			if (i >= argc) usage();
			dr[2] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-soft") == 0) {
			i++;
			if (i >= argc) usage();
			toplevelsoftening = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-mass") == 0) {
			i++;
			if (i >= argc) usage();
			toplevelmass = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-LUsf") == 0) {
			i++;
			if (i >= argc) usage();
			LUsf = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-VUsf") == 0) {
			i++;
			if (i >= argc) usage();
			VUsf = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-MUsf") == 0) {
			i++;
			if (i >= argc) usage();
			MUsf = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-b2dmfac") == 0) {
			i++;
			if (i >= argc) usage();
			b2dmscalefac = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-noscaling") == 0) {
			dr[0] = 0;
			dr[1] = 0;
			dr[2] = 0;
			LUsf = 1;
			VUsf = 1;
			MUsf = 1;
			b2dmscalefac = 1;
			toplevelmass = -1;
			toplevelsoftening = -2;
			i++;
			}
		else if (strcmp(argv[i],"-softfac") == 0) {
			i++;
			if (i >= argc) usage();
			softfac = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-refstep") == 0) {
			i++;
			if (i >= argc) usage();
			refinementstep = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-pos") == 0) {
			i++;
			if (i >= argc) usage();
			strcpy(PosFileName,argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-vel") == 0) {
			i++;
			if (i >= argc) usage();
			strcpy(VelFileName,argv[i]);
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

	if (toplevelmass == -1.0) mrmassfromfile = 1;

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
	N1Dlow = PosHeader.dims[0];
	dx = PosManifest.dx;
	Lmax = PosHeader.Lmax;
	LBox = N1Dlow*dx;

	acurrent = PosHeader.aBegin;
	h100 = PosManifest.h100;
	OmegaM0 = PosManifest.OmegaX+PosManifest.OmegaB;
	OmegaDM0 = PosManifest.OmegaX;
	OmegaB0 = PosManifest.OmegaB;
	OmegaL0 = PosManifest.OmegaL;
	OmegaK0 = 0;
	OmegaR0 = PosManifest.OmegaN;

	th.time = acurrent;
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

	TU_TIPSY_DEFAULT = (H_TIPSY_DEFAULT*1000)/(h100*100); /* kpc km^-1 s */
	LU_TIPSY_DEFAULT = LBox*1000/h100; /* kpc */
	VU_TIPSY_DEFAULT = LU_TIPSY_DEFAULT/TU_TIPSY_DEFAULT; /* km s^-1 */
	MU_TIPSY_DEFAULT = rhocrit*h100*h100*LU_TIPSY_DEFAULT*LU_TIPSY_DEFAULT*LU_TIPSY_DEFAULT; /* Mo */

	LU_GIC = 1000/h100; /* kpc */
	VU_GIC = 1; /* km s^-1 */
	MU_GIC = 1; /* Mo */
	TU_GIC = LU_GIC/VU_GIC; /* kpc km^-1 s */

	if (LUsf < 0) LUsf = 1.0/LBox; /* LU_GIC / LU_TIPSY_DEFAULT */
	if (VUsf < 0) VUsf = VU_GIC/VU_TIPSY_DEFAULT;
	if (MUsf < 0) MUsf = MU_GIC/MU_TIPSY_DEFAULT;
	if (b2dmscalefac < 0) b2dmscalefac = OmegaM0/OmegaDM0;
	csvelscalefac = 1/acurrent;

	LU_TIPSY = LU_GIC/LUsf; /* kpc */
	VU_TIPSY = VU_GIC/VUsf; /* km s^-1 */
	MU_TIPSY = MU_GIC/MUsf; /* Mo */
	TU_TIPSY = LU_TIPSY/VU_TIPSY; /* kpc km^-1 s */

	/*
	** Get output file ready
	*/

	xdrstdio_create(&xdrs,stdout,XDR_ENCODE);
	write_tipsy_xdr_header(&xdrs,&th);

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

		if (toplevelsoftening < 0) toplevelsoftening = (LBox/(N1Dlow*softfac))*(1000/(h100*LU_TIPSY));
		if (toplevelmass < 0) toplevelmass = OmegaM0/Ntot;
		soft = toplevelsoftening;
		mass = toplevelmass;

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
						tgp.pos[0] = PosXRec[j]*LUsf + dr[0];
						tgp.pos[1] = PosYRec[j]*LUsf + dr[1];
						tgp.pos[2] = PosZRec[j]*LUsf + dr[2];
						tgp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
						tgp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
						tgp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
						tgp.mass = mass;
						tgp.rho = 0;
						tgp.temp = 0;
						tgp.hsmooth = soft;
						tgp.metals = 0;
						tgp.phi = 0;
						if (index < Ntot) write_tipsy_xdr_gas(&xdrs,&tgp);
						}
					else if (particletype == 1) {
						tdp.pos[0] = PosXRec[j]*LUsf + dr[0];
						tdp.pos[1] = PosYRec[j]*LUsf + dr[1];
						tdp.pos[2] = PosZRec[j]*LUsf + dr[2];
						tdp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
						tdp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
						tdp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
						tdp.mass = mass;
						tdp.eps = soft;
						tdp.phi = 0;
						if (index < Ntot) write_tipsy_xdr_dark(&xdrs,&tdp);
						}
					else {
						tsp.pos[0] = PosXRec[j]*LUsf + dr[0];
						tsp.pos[1] = PosYRec[j]*LUsf + dr[1];
						tsp.pos[2] = PosZRec[j]*LUsf + dr[2];
						tsp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
						tsp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
						tsp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
						tsp.mass = mass;
						tsp.metals = 0;
						tsp.tform = 0;
						tsp.eps = soft;
						tsp.phi = 0;
						if (index < Ntot) write_tipsy_xdr_star(&xdrs,&tsp);
						}
					}
				else {
					if (particletype == 0) {
						tgpdpp.pos[0] = PosXRec[j]*LUsf + dr[0];
						tgpdpp.pos[1] = PosYRec[j]*LUsf + dr[1];
						tgpdpp.pos[2] = PosZRec[j]*LUsf + dr[2];
						tgpdpp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
						tgpdpp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
						tgpdpp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
						tgpdpp.mass = mass;
						tgpdpp.rho = 0;
						tgpdpp.temp = 0;
						tgpdpp.hsmooth = soft;
						tgpdpp.metals = 0;
						tgpdpp.phi = 0;
						if (index < Ntot) write_tipsy_xdr_gas_dpp(&xdrs,&tgpdpp);
						}
					else if (particletype == 1) {
						tdpdpp.pos[0] = PosXRec[j]*LUsf + dr[0];
						tdpdpp.pos[1] = PosYRec[j]*LUsf + dr[1];
						tdpdpp.pos[2] = PosZRec[j]*LUsf + dr[2];
						tdpdpp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
						tdpdpp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
						tdpdpp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
						tdpdpp.mass = mass;
						tdpdpp.eps = soft;
						tdpdpp.phi = 0;
						if (index < Ntot) write_tipsy_xdr_dark_dpp(&xdrs,&tdpdpp);
						}
					else {
						tspdpp.pos[0] = PosXRec[j]*LUsf + dr[0];
						tspdpp.pos[1] = PosYRec[j]*LUsf + dr[1];
						tspdpp.pos[2] = PosZRec[j]*LUsf + dr[2];
						tspdpp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
						tspdpp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
						tspdpp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
						tspdpp.mass = mass;
						tspdpp.metals = 0;
						tspdpp.tform = 0;
						tspdpp.eps = soft;
						tspdpp.phi = 0;
						if (index < Ntot) write_tipsy_xdr_star_dpp(&xdrs,&tspdpp);
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

		if (verboselevel > 0) {
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
				if (toplevelsoftening < 0) toplevelsoftening = (LBox/(N1Dlow*softfac))*(1000/(h100*LU_TIPSY));
				soft = toplevelsoftening;
				}
			else {
				soft = toplevelsoftening/refinementstep;
				}
			if (k == 0) { /* top level mass */
				if (toplevelmass < 0) toplevelmass = OmegaM0/Nlow;
				mass = toplevelmass;
				}
			else {
				mass = toplevelmass/pow(refinementstep,3.0);
				}
			if (mrmassfromfile == 1) {
				assert(PosLevelHeader[k].Mlev == VelLevelHeader[k].Mlev);
				mass = PosLevelHeader[k].Mlev*b2dmscalefac*MUsf;
				if (k == 0) toplevelmass = mass;
				}

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
							tgp.pos[0] = PosXRec[j]*LUsf + dr[0];
							tgp.pos[1] = PosYRec[j]*LUsf + dr[1];
							tgp.pos[2] = PosZRec[j]*LUsf + dr[2];
							tgp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
							tgp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
							tgp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
							tgp.mass = mass;
							tgp.rho = 0;
							tgp.temp = 0;
							tgp.hsmooth = soft;
							tgp.metals = 0;
							tgp.phi = 0;
							if (index < Nlev) write_tipsy_xdr_gas(&xdrs,&tgp);
							}
						else if (particletype == 1) {
							tdp.pos[0] = PosXRec[j]*LUsf + dr[0];
							tdp.pos[1] = PosYRec[j]*LUsf + dr[1];
							tdp.pos[2] = PosZRec[j]*LUsf + dr[2];
							tdp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
							tdp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
							tdp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
							tdp.mass = mass;
							tdp.eps = soft;
							tdp.phi = 0;
							if (index < Nlev) write_tipsy_xdr_dark(&xdrs,&tdp);
							}
						else {
							tsp.pos[0] = PosXRec[j]*LUsf + dr[0];
							tsp.pos[1] = PosYRec[j]*LUsf + dr[1];
							tsp.pos[2] = PosZRec[j]*LUsf + dr[2];
							tsp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
							tsp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
							tsp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
							tsp.mass = mass;
							tsp.metals = 0;
							tsp.tform = 0;
							tsp.eps = soft;
							tsp.phi = 0;
							if (index < Nlev) write_tipsy_xdr_star(&xdrs,&tsp);
							}
						}
					else {
						if (particletype == 0) {
							tgpdpp.pos[0] = PosXRec[j]*LUsf + dr[0];
							tgpdpp.pos[1] = PosYRec[j]*LUsf + dr[1];
							tgpdpp.pos[2] = PosZRec[j]*LUsf + dr[2];
							tgpdpp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
							tgpdpp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
							tgpdpp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
							tgpdpp.mass = mass;
							tgpdpp.rho = 0;
							tgpdpp.temp = 0;
							tgpdpp.hsmooth = soft;
							tgpdpp.metals = 0;
							tgpdpp.phi = 0;
							if (index < Nlev) write_tipsy_xdr_gas_dpp(&xdrs,&tgpdpp);
							}
						else if (particletype == 1) {
							tdpdpp.pos[0] = PosXRec[j]*LUsf + dr[0];
							tdpdpp.pos[1] = PosYRec[j]*LUsf + dr[1];
							tdpdpp.pos[2] = PosZRec[j]*LUsf + dr[2];
							tdpdpp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
							tdpdpp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
							tdpdpp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
							tdpdpp.mass = mass;
							tdpdpp.eps = soft;
							tdpdpp.phi = 0;
							if (index < Nlev) write_tipsy_xdr_dark_dpp(&xdrs,&tdpdpp);
							}
						else {
							tspdpp.pos[0] = PosXRec[j]*LUsf + dr[0];
							tspdpp.pos[1] = PosYRec[j]*LUsf + dr[1];
							tspdpp.pos[2] = PosZRec[j]*LUsf + dr[2];
							tspdpp.vel[0] = VelXRec[j]*csvelscalefac*VUsf;
							tspdpp.vel[1] = VelYRec[j]*csvelscalefac*VUsf;
							tspdpp.vel[2] = VelZRec[j]*csvelscalefac*VUsf;
							tspdpp.mass = mass;
							tspdpp.metals = 0;
							tspdpp.tform = 0;
							tspdpp.eps = soft;
							tspdpp.phi = 0;
							if (index < Nlev) write_tipsy_xdr_star_dpp(&xdrs,&tspdpp);
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

			if (verboselevel > 0) {
				fprintf(stderr,"L %d Lmax %d Nlev %ld Softening %.6e LU_TIPSY = %.6e kpc Mass %.6e MU_TIPSY = %.6e Mo\n",PosLevelHeader[k].L,PosLevelHeader[k].Lmax,Nlev,soft,soft*LU_TIPSY,mass,mass*MU_TIPSY);
				if (k == Lmax) fprintf(stderr,"\n");
				}
			}
		}

	xdr_destroy(&xdrs);

	/*
	** Write out some additional stuff depending on verbose level
	*/

	if (verboselevel > 0) {
		fprintf(stderr,"Parameters from GIC file:\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"Name    : %s\n",PosManifest.name);
		fprintf(stderr,"OmegaB  : %.6e\n",PosManifest.OmegaB);
		fprintf(stderr,"OmegaDM : %.6e\n",PosManifest.OmegaX);
		fprintf(stderr,"OmegaL  : %.6e\n",PosManifest.OmegaL);
		fprintf(stderr,"OmegaN  : %.6e\n",PosManifest.OmegaN);
		fprintf(stderr,"h100    : %.6e\n",h100);
		fprintf(stderr,"Deltax  : %.6e chimp\n",dx);
		fprintf(stderr,"ns      : %.6e\n",PosManifest.ns);
		fprintf(stderr,"sigma8  : %.6e\n",PosManifest.s8);
		fprintf(stderr,"kp      : %.6e\n",PosManifest.kp);
		fprintf(stderr,"aBegin  : %.6e\n",PosHeader.aBegin);
		fprintf(stderr,"DeltaDC : %.6e\n",PosHeader.DeltaDC);
		fprintf(stderr,"NX      : %d\n",PosHeader.dims[0]);
		fprintf(stderr,"NY      : %d\n",PosHeader.dims[1]);
		fprintf(stderr,"NZ      : %d\n",PosHeader.dims[2]);
		fprintf(stderr,"Seed    : %d\n",PosHeader.seed);
		fprintf(stderr,"Nrec    : %d\n",Nrec);
		fprintf(stderr,"LBox    : %.6e chimp\n",LBox);
		fprintf(stderr,"Lmax    : %d\n",Lmax);
		fprintf(stderr,"\n");
		fprintf(stderr,"Cosmology:\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"OmegaM0  : %.6e\n",OmegaM0);
		fprintf(stderr,"OmegaDM0 : %.6e\n",OmegaDM0);
		fprintf(stderr,"OmegaB0  : %.6e\n",OmegaB0);
		fprintf(stderr,"OmegaL0  : %.6e\n",OmegaL0);
		fprintf(stderr,"OmegaK0  : %.6e\n",OmegaK0);
		fprintf(stderr,"OmegaR0  : %.6e\n",OmegaR0);
		fprintf(stderr,"h100     : %.6e\n",h100);
		fprintf(stderr,"\n");
		fprintf(stderr,"Used values:\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"drx       : %.6e LU_TIPSY\n",dr[0]);
		fprintf(stderr,"dry       : %.6e LU_TIPSY\n",dr[1]);
		fprintf(stderr,"drz       : %.6e LU_TIPSY\n",dr[2]);
		fprintf(stderr,"LUsf      : %.6e\n",LUsf);
		fprintf(stderr,"VUsf      : %.6e\n",VUsf);
		fprintf(stderr,"MUsf      : %.6e\n",MUsf);
		fprintf(stderr,"softfac   : %.6e\n",softfac);
		fprintf(stderr,"Softening : %.6e LU_TIPSY (toplevel)\n",toplevelsoftening);
		fprintf(stderr,"Mass      : %.6e MU_TIPSY (toplevel)\n",toplevelmass);
		fprintf(stderr,"\n");
		fprintf(stderr,"Resulting internal GIC units:\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"LU_GIC : %.6e kpc\n",LU_GIC);
		fprintf(stderr,"TU_GIC : %.6e Gyr\n",TU_GIC/VelConvertFac);
		fprintf(stderr,"VU_GIC : %.6e km s^-1 = %.6e kpc Gyr^-1\n",VU_GIC,VU_GIC*VelConvertFac);
		fprintf(stderr,"MU_GIC : %.6e Mo\n",MU_GIC);
		fprintf(stderr,"\n");
		fprintf(stderr,"Resulting internal tipsy units:\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"LU_TIPSY : %.6e kpc\n",LU_TIPSY);
		fprintf(stderr,"TU_TIPSY : %.6e Gyr\n",TU_TIPSY/VelConvertFac);
		fprintf(stderr,"VU_TIPSY : %.6e km s^-1 = %.6e kpc Gyr^-1\n",VU_TIPSY,VU_TIPSY*VelConvertFac);
		fprintf(stderr,"MU_TIPSY : %.6e Mo\n",MU_TIPSY);
		fprintf(stderr,"\n");
		fprintf(stderr,"Time: %g Ntotal: %u Ngas: %u Ndark: %u Nstar: %u\n",
			th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
		}

	exit(0);
	}

void usage(void) {

	fprintf(stderr,"\n");
	fprintf(stderr,"%s (%s)\n",NAME,VERSION);
	fprintf(stderr,"\n");
	fprintf(stderr,"Program converts GIC native binary format to tipsy XDR format.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Please specify the following parameters:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-spp             : set this flag if output file has single precision positions (default)\n");
	fprintf(stderr,"-dpp             : set this flag if output file has double precision positions\n");
	fprintf(stderr,"-sr              : set this flag if input only has single resolution particles (default)\n");
	fprintf(stderr,"-mr              : set this flag if input has multi resolution particles\n");
	fprintf(stderr,"-gas             : set this flag if input particles are gas particles\n");
	fprintf(stderr,"-dark            : set this flag if input particles are dark particles (default)\n");
	fprintf(stderr,"-star            : set this flag if input particles are star particles\n");
	fprintf(stderr,"-drx <value>     : shift along x-axis [LU_TIPSY] (default: -0.5 LU_TIPSY)\n");
	fprintf(stderr,"-dry <value>     : shift along y-axis [LU_TIPSY] (default: -0.5 LU_TIPSY)\n");
	fprintf(stderr,"-drz <value>     : shift along z-axis [LU_TIPSY] (default: -0.5 LU_TIPSY)\n");
	fprintf(stderr,"-soft <value>    : softening length of top level particles [LU_TIPSY] (default: 1/softfac mean particle separation => LBox/(N1Dlow*softfac) LU_TIPSY)\n");
	fprintf(stderr,"-mass <value>    : mass of top level particles [MU_TIPSY] (default: OmegaM0/Nlow - if you want masses from file in mr case set -1)\n");
	fprintf(stderr,"-LUsf <value>    : length unit scale factor (default: LU_GIC/LU_TIPSY_DEFAULT)\n");
	fprintf(stderr,"-VUsf <value>    : velocity unit scale factor (default: VU_GIC/VU_TIPSY_DEFAULT)\n");
	fprintf(stderr,"-MUsf <value>    : mass unit scale factor (default: MU_GIC/MU_TIPSY_DEFAULT - only in mr case when read from file)\n");
	fprintf(stderr,"-b2dmfac <value> : baryon to dark matter scale factor (default: OmegaM0/OmegaDM0)\n");
	fprintf(stderr,"-softfac <value> : softening factor (default: 50)\n");
	fprintf(stderr,"-refstep <value> : refinement step factor (default: 2)\n");
	fprintf(stderr,"-noscaling       : no scaling of data\n");
	fprintf(stderr,"-verbose         : verbose\n");
	fprintf(stderr,"-pos <name>      : positions input file GIC native binary format\n");
	fprintf(stderr,"-vel <name>      : velocities input file GIC native binary format\n");
	fprintf(stderr,"> <name>         : output file in tipsy XDR format\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Other options:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-h or -help : display this help and exit\n");
	fprintf(stderr,"-version    : display version information and exit\n");
	fprintf(stderr,"\n");
	exit(1);
	}

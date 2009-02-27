/* 
** artb2ts.c
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

void usage(void);

int main(int argc, char **argv) {

    int i = 0, j = 0, k = 0;
    int header, trailer;
    int positionprecision, verboselevel, doswap, mrmassfromfile;
    long Ntot, index, N[10];
    int Nrec, Npage, N1Dlow, Nlow, L, Lmax;
    double posscalefac, velscalefac, massscalefac, refinementstep;
    double particlemass, particlesoftening;
    double toplevelmass, toplevelsoftening;
    double dr[3], aBegin, VelConvertFac, rhocrit, LBox, softfac;
    double OmegaM0, OmegaDM0, OmegaB0, OmegaL0, OmegaK0, OmegaR0, h100;
    double H_Tipsy, TU_Tipsy, LU_Tipsy, VU_Tipsy, MU_Tipsy;
    double TU_ART, LU_ART, VU_ART, MU_ART;
    double mass[10], soft[10];
    char HeaderFileName[256], DataFileName[256];
    char banner[45];
    float rx, ry, rz, vx, vy, vz;
    TIPSY_HEADER th;
    DARK_PARTICLE dp;
    DARK_PARTICLE_DPP dpdpp;
    ART_HEADER ah;
    FILE *HeaderFile;
    FILE *PosXFile, *PosYFile, *PosZFile;
    FILE *VelXFile, *VelYFile, *VelZFile;
    XDR xdrs;

    /*
    ** Set some default values
    */

    positionprecision = 0;
    verboselevel = 0;
    doswap = 0;
    mrmassfromfile = 0;
    dr[0] = -0.5;
    dr[1] = -0.5;
    dr[2] = -0.5;
    refinementstep = 2;
    particlesoftening = -2;
    particlemass = -2;
    posscalefac = -2;
    velscalefac = -2;
    massscalefac = -2;
    softfac = 20;
    toplevelsoftening = 0;
    toplevelmass = 0;
    Nrec = 1024*1024;
    H_Tipsy = sqrt(8*M_PI/3); /* TU_Tipsy^-1 */
    VelConvertFac = 1.0227121651152353693;
    rhocrit = 277.53662719; /* h^2 Mo kpc^-3 */
    for (i = 0; i < 10; i++) {
	N[i] = 0;
	mass[i] = 0;
	soft[i] = 0;
	}

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
        else if (strcmp(argv[i],"-massfac") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            massscalefac = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-softfac") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            softfac = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-refstep") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            refinementstep = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-header") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            strcpy(HeaderFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-data") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            strcpy(DataFileName,argv[i]);
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

    if (particlemass == -1.0) mrmassfromfile = 1;

    /*
    ** Read heder file
    */

    HeaderFile = fopen(HeaderFileName,"r");
    assert(HeaderFile != NULL);
    assert(fread(&header,sizeof(int),1,HeaderFile) == 1);
    if (header != 45+sizeof(ART_HEADER)) {
	doswap = 1;
	flip_4byte(&header,sizeof(int),1);
	}
    assert(header == 45+sizeof(ART_HEADER));
    assert(fread(&banner,sizeof(char),45,HeaderFile) == 45);
    assert(fread(&ah,sizeof(ART_HEADER),1,HeaderFile) == 1);
    if (doswap) flip_4byte(&ah,sizeof(ART_HEADER),1);
    assert(fread(&trailer,sizeof(int),1,HeaderFile) == 1);
    if (doswap) flip_4byte(&trailer,sizeof(int),1);
    fclose(HeaderFile);

    LBox = ah.fill[79];
    Lmax = ah.Nspecies-1;
    Ntot = ah.num[Lmax];
    Npage = (Ntot+Nrec-1)/Nrec;
    N1Dlow = ah.Ngrid;
    Nlow = N1Dlow*N1Dlow*N1Dlow;

    aBegin = ah.aexpn;
    h100 = ah.hubble;
    OmegaM0 = ah.Om0;
    OmegaDM0 = ah.Om0-ah.Omb0;
    OmegaB0 = ah.Omb0;
    OmegaL0 = ah.Oml0;
    OmegaK0 = ah.Ocurv;
    OmegaR0 = 0;

    /*
    ** Calculate scaling factors
    */

    TU_Tipsy = (H_Tipsy*1000)/(h100*100); /* kpc km^-1 s */
    LU_Tipsy = LBox*1000/h100; /* kpc */
    VU_Tipsy = LU_Tipsy/TU_Tipsy; /* km s^-1 */
    MU_Tipsy = rhocrit*h100*h100*LU_Tipsy*LU_Tipsy*LU_Tipsy; /* Mo */

    TU_ART = 1000*2.0/(100*h100*sqrt(OmegaM0)); /* kpc km^-1 s */
    LU_ART = LBox*1000/(h100*ah.Ngrid); /* kpc */
    VU_ART = LU_ART/TU_ART; /* km s^-1 */
    MU_ART = rhocrit*h100*h100*LU_ART*LU_ART*LU_ART*OmegaM0; /* Mo */

    if (posscalefac < 0) posscalefac = 1.0/ah.Ngrid; /* LU_ART / LU_Tipsy */
    if (velscalefac < 0) velscalefac = VU_ART/(aBegin*VU_Tipsy);
    if (massscalefac < 0) massscalefac = MU_ART/MU_Tipsy;
    
    if (particlesoftening < 0) particlesoftening = 1.0/(N1Dlow*softfac);
    toplevelsoftening = particlesoftening;
    if (particlemass < 0) particlemass = OmegaM0/Nlow;
    if (mrmassfromfile == 1) particlemass = ah.mass[Lmax]*(OmegaM0/OmegaDM0)*massscalefac;
    toplevelmass = particlemass;

    /*
    ** Get data files ready
    */

    PosXFile = fopen(DataFileName,"r");
    assert(PosXFile != NULL);
    PosYFile = fopen(DataFileName,"r");
    assert(PosYFile != NULL);
    PosZFile = fopen(DataFileName,"r");
    assert(PosZFile != NULL);
    VelXFile = fopen(DataFileName,"r");
    assert(VelXFile != NULL);
    VelYFile = fopen(DataFileName,"r");
    assert(VelYFile != NULL);
    VelZFile = fopen(DataFileName,"r");
    assert(VelZFile != NULL);

    /*
    ** Get output file ready
    */

    th.time = aBegin;
    th.ntotal = Ntot;
    th.ndim = 3;
    th.ngas = 0; 
    th.ndark = Ntot;
    th.nstar = 0;

    xdrstdio_create(&xdrs,stdout,XDR_ENCODE);
    write_tipsy_standard_header(&xdrs,&th);

    /*
    ** Read an process data
    */

    index = 0;
    for (i = 0; i < Npage; i++) {

	/*
	** Get file pointers ready
	*/

	assert(fseek(PosYFile,1*Nrec*sizeof(float),SEEK_CUR) == 0);
	assert(fseek(PosZFile,2*Nrec*sizeof(float),SEEK_CUR) == 0);
	assert(fseek(VelXFile,3*Nrec*sizeof(float),SEEK_CUR) == 0);
	assert(fseek(VelYFile,4*Nrec*sizeof(float),SEEK_CUR) == 0);
	assert(fseek(VelZFile,5*Nrec*sizeof(float),SEEK_CUR) == 0);

	/*
	** Go through records
	*/

	for (j = 0; j < Nrec; j++) {

	    assert(fread(&rx,sizeof(float),1,PosXFile) == 1);
	    assert(fread(&ry,sizeof(float),1,PosYFile) == 1);
	    assert(fread(&rz,sizeof(float),1,PosZFile) == 1);
	    assert(fread(&vx,sizeof(float),1,VelXFile) == 1);
	    assert(fread(&vy,sizeof(float),1,VelYFile) == 1);
	    assert(fread(&vz,sizeof(float),1,VelZFile) == 1);

	    if (doswap) flip_4byte(&rx,sizeof(float),1);
	    if (doswap) flip_4byte(&ry,sizeof(float),1);
	    if (doswap) flip_4byte(&rz,sizeof(float),1);
	    if (doswap) flip_4byte(&vx,sizeof(float),1);
	    if (doswap) flip_4byte(&vy,sizeof(float),1);
	    if (doswap) flip_4byte(&vz,sizeof(float),1);

            /*
	    ** Set particle properties
            */

	    for (k = Lmax; k >=0; k--) {
		L = Lmax-k;
		if (ah.num[k] > index) {
		    particlesoftening = toplevelsoftening/pow(refinementstep,L);
		    particlemass = toplevelmass/pow(refinementstep,3.0*L);
		    if (mrmassfromfile == 1) particlemass = ah.mass[k]*(OmegaM0/OmegaDM0)*massscalefac;
		    if (N[L] == 0) N[L] = ah.num[k];
		    if (mass[L] == 0) mass[L] = particlemass;
		    if (soft[L] == 0) soft[L] = particlesoftening;
		    }
		}
	    if (positionprecision == 0) {
		dp.pos[0] = (rx-1.0)*posscalefac + dr[0];
		dp.pos[1] = (ry-1.0)*posscalefac + dr[1];
		dp.pos[2] = (rz-1.0)*posscalefac + dr[2];
		dp.vel[0] = vx*velscalefac;
		dp.vel[1] = vy*velscalefac;
		dp.vel[2] = vz*velscalefac;
		dp.mass = particlemass;
		dp.eps = particlesoftening;
		dp.phi = 0;
		if (index < Ntot) write_tipsy_standard_dark(&xdrs,&dp);
		}
	    else {
		dpdpp.pos[0] = (rx-1.0)*posscalefac + dr[0];
		dpdpp.pos[1] = (ry-1.0)*posscalefac + dr[1];
		dpdpp.pos[2] = (rz-1.0)*posscalefac + dr[2];
		dpdpp.vel[0] = vx*velscalefac;
		dpdpp.vel[1] = vy*velscalefac;
		dpdpp.vel[2] = vz*velscalefac;
		dpdpp.mass = particlemass;
		dpdpp.eps = particlesoftening;
		dpdpp.phi = 0;
		if (index < Ntot) write_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
		}
	    index++;
	    }

	/*
	** Move all pointers to the end of the record
	*/

	assert(fseek(PosXFile,5*Nrec*sizeof(float),SEEK_CUR) == 0);
	assert(fseek(PosYFile,4*Nrec*sizeof(float),SEEK_CUR) == 0);
	assert(fseek(PosZFile,3*Nrec*sizeof(float),SEEK_CUR) == 0);
	assert(fseek(VelXFile,2*Nrec*sizeof(float),SEEK_CUR) == 0);
	assert(fseek(VelYFile,1*Nrec*sizeof(float),SEEK_CUR) == 0);
	}

    xdr_destroy(&xdrs);

    /*
    ** Write out some additional stuff depending on verbose level
    */

    if (verboselevel >= 1) {
	fprintf(stderr,"There are %d refinement levels:\n\n",Lmax+1);
	for (L = 0; L <= Lmax; L++) {
	    fprintf(stderr,"L %d Lmax %d Nlev %ld Softening %.6e LU = %.6e kpc Mass %.6e MU = %.6e Mo\n",L,Lmax,N[L],soft[L],soft[L]*LU_Tipsy,mass[L],mass[L]*MU_Tipsy);
	    if (L == Lmax) fprintf(stderr,"\n");
	    }
	fprintf(stderr,"Parameters from ART file:\n\n");
        fprintf(stderr,"aexpn    : %.6e\n",ah.aexpn);
        fprintf(stderr,"aexp0    : %.6e\n",ah.aexp0);
        fprintf(stderr,"amplt    : %.6e\n",ah.amplt);
        fprintf(stderr,"astep    : %.6e\n",ah.astep);
        fprintf(stderr,"istep    : %d\n",ah.istep);
        fprintf(stderr,"partw    : %.6e\n",ah.partw);
        fprintf(stderr,"tintg    : %.6e\n",ah.tintg);
        fprintf(stderr,"ekin     : %.6e\n",ah.ekin);
        fprintf(stderr,"ekin1    : %.6e\n",ah.ekin1);
        fprintf(stderr,"ekin2    : %.6e\n",ah.ekin2);
        fprintf(stderr,"au0      : %.6e\n",ah.au0);
        fprintf(stderr,"aeu0     : %.6e\n",ah.aeu0);
        fprintf(stderr,"Nrow     : %d\n",ah.Nrow);
        fprintf(stderr,"Ngrid    : %d\n",ah.Ngrid);
        fprintf(stderr,"Nspecies : %d\n",ah.Nspecies);
        fprintf(stderr,"Nseed    : %d\n",ah.Nseed);
        fprintf(stderr,"Om0      : %.6e\n",ah.Om0);
        fprintf(stderr,"Oml0     : %.6e\n",ah.Oml0);
        fprintf(stderr,"hubble   : %.6e\n",ah.hubble);
        fprintf(stderr,"Wp5      : %.6e\n",ah.Wp5);
        fprintf(stderr,"Ocurv    : %.6e\n",ah.Ocurv);
        fprintf(stderr,"Omb0     : %.6e\n",ah.Omb0);
        fprintf(stderr,"Banner   : %s\n",banner);
        fprintf(stderr,"Nrec     : %d\n",Nrec);
        fprintf(stderr,"LBox     : %.6e chimp\n",LBox);
        fprintf(stderr,"Lmax     : %d\n\n",Lmax);
        fprintf(stderr,"Cosmology:\n\n");
        fprintf(stderr,"OmegaM0  : %.6e\n",OmegaM0);
        fprintf(stderr,"OmegaDM0 : %.6e\n",OmegaDM0);
        fprintf(stderr,"OmegaB0  : %.6e\n",OmegaB0);
        fprintf(stderr,"OmegaL0  : %.6e\n",OmegaL0);
        fprintf(stderr,"OmegaK0  : %.6e\n",OmegaK0);
        fprintf(stderr,"OmegaR0  : %.6e\n",OmegaR0);
        fprintf(stderr,"h100     : %.6e\n\n",h100);
	fprintf(stderr,"Used values:\n\n");
        fprintf(stderr,"drx       : %.6e LU\n",dr[0]);
        fprintf(stderr,"dry       : %.6e LU\n",dr[1]);
        fprintf(stderr,"drz       : %.6e LU\n",dr[2]);
        fprintf(stderr,"posfac    : %.6e\n",posscalefac);
        fprintf(stderr,"velfac    : %.6e\n",velscalefac);
        fprintf(stderr,"massfac   : %.6e\n",massscalefac);
        fprintf(stderr,"softfac   : %.6e\n",softfac);
        fprintf(stderr,"Softening : %.6e LU (toplevel)\n",toplevelsoftening);
        fprintf(stderr,"Mass      : %.6e MU (toplevel)\n\n",toplevelmass);
        fprintf(stderr,"Resulting internal ART units:\n\n");
        fprintf(stderr,"LU : %.6e kpc\n",LU_ART);
        fprintf(stderr,"TU : %.6e Gyr\n",TU_ART/VelConvertFac);
        fprintf(stderr,"VU : %.6e km s^-1 = %.6e kpc Gyr^-1\n",VU_ART,VU_ART*VelConvertFac);
        fprintf(stderr,"MU : %.6e Mo\n\n",MU_ART);
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
    fprintf(stderr,"Program converts ART binary format to tipsy standard binary format.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp             : set this flag if output file has single precision positions (default)\n");
    fprintf(stderr,"-dpp             : set this flag if output file has double precision positions\n");
    fprintf(stderr,"-drx <value>     : shift along x-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-dry <value>     : shift along y-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-drz <value>     : shift along z-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-soft <value>    : softening length of top level particles [LU] (default: 1/softfac mean particle separation => 1/[Nlow^{-3}*softfac] LU)\n");
    fprintf(stderr,"-mass <value>    : mass of top level particles [MU] (default: OmegaM0/Nlow - if you want masses from file in mr case set -1)\n");
    fprintf(stderr,"-posfac <value>  : position scale factor (default: LU_ART/LU_Tipsy)\n");
    fprintf(stderr,"-velfac <value>  : velocity scale factor (default: VU_ART/[a*VU_Tipsy] where a is the scale factor)\n");
    fprintf(stderr,"-massfac <value> : mass scale factor (default: MU_ART/MU_Tipsy)\n");
    fprintf(stderr,"-softfac <value> : softening factor (default: 20)\n");
    fprintf(stderr,"-refstep <value> : refinement step factor (default: 2)\n");
    fprintf(stderr,"-v               : more informative output to screen\n");
    fprintf(stderr,"-header <name>   : header input file in ART binary format\n");
    fprintf(stderr,"-data <name>     : data input file in ART binary format\n");
    fprintf(stderr,"> <name>         : output file in tipsy standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }
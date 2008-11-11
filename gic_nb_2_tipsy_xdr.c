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
#include <iof.h>

void usage(void);

int main(int argc, char **argv) {

    int i, j, k, l;
    int positionprecision, verboselevel;
    int header,trailer;
    int NX, NY, NZ;
    float OmB, OmX, OmL, OmN, h100, cell, enn, sigma8, akpivot, aBegin, delDC;
    float tempfloat;
    double posscalefac, velscalefac, particlemass, particlesoftening;
    double dr[3], H_Tipsy, TU_Tipsy, LU_Tipsy, VU_Tipsy, MU_Tipsy, VelConvertFac, rhocrit;
    char name[256];
    TIPSY_HEADER th;
    TIPSY_STRUCTURE *ts;
    TIPSY_STRUCTURE_DPP *tsdpp;

    ts = NULL;
    tsdpp = NULL;
    positionprecision = 0;
    verboselevel = 0;
    dr[0] = -0.5;
    dr[1] = -0.5;
    dr[2] = -0.5;
    particlesoftening = 0;
    H_Tipsy = sqrt(8*M_PI/3); /* TU_Tipsy^-1 */
    VelConvertFac = 1.0227121651152353693;
    rhocrit = 277.53662719; /* h^2 Mo kpc^-3 */
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
    ** Read in raw data first
    */

    assert(fread(&header,sizeof(int),1,stdin) == 1);
    assert(fread(name,sizeof(name),1,stdin) == 1);
    assert(fread(&OmB,sizeof(float),1,stdin) == 1);
    assert(fread(&OmX,sizeof(float),1,stdin) == 1);
    assert(fread(&OmL,sizeof(float),1,stdin) == 1);
    assert(fread(&OmN,sizeof(float),1,stdin) == 1);
    assert(fread(&h100,sizeof(float),1,stdin) == 1);
    assert(fread(&cell,sizeof(float),1,stdin) == 1);
    assert(fread(&enn,sizeof(float),1,stdin) == 1);
    assert(fread(&sigma8,sizeof(float),1,stdin) == 1);
    assert(fread(&akpivot,sizeof(float),1,stdin) == 1);
    assert(fread(&trailer,sizeof(int),1,stdin) == 1);
    assert(header == trailer);

    assert(fread(&header,sizeof(int),1,stdin) == 1);
    assert(fread(&aBegin,sizeof(float),1,stdin) == 1);
    assert(fread(&delDC,sizeof(float),1,stdin) == 1);
    assert(fread(&NX,sizeof(int),1,stdin) == 1);
    assert(fread(&NY,sizeof(int),1,stdin) == 1);
    assert(fread(&NZ,sizeof(int),1,stdin) == 1);
    assert(fread(&trailer,sizeof(int),1,stdin) == 1);
    assert(header == trailer);

    th.time = aBegin;
    th.ntotal = NX*NY*NZ;
    th.ndim = 3;
    th.ngas = 0;
    th.ndark = NX*NY*NZ;
    th.nstar = 0;

    if (positionprecision == 0) {
	tsdpp = NULL;
	ts = malloc(sizeof(TIPSY_STRUCTURE));
	assert(ts != NULL);
	ts->th = &th;
	ts->gp = NULL;
	ts->dp = malloc(NX*NY*NZ*sizeof(DARK_PARTICLE));
	assert(ts->dp != NULL);
	ts->sp = NULL;
	for (i = 0; i < 3; i++) {
	    l = 0;
	    for(j = 0; j < NZ; j++) {
		assert(fread(&header,sizeof(int),1,stdin) == 1);
		for (k = 0; k < NX*NY; k++) {
		    assert(fread(&tempfloat,sizeof(float),1,stdin) == 1);
		    ts->dp[l].vel[i] = tempfloat;
		    l++;
		    }
		assert(fread(&trailer,sizeof(int),1,stdin) == 1);
		assert(header == trailer);
		}
	    l = 0;
            for(j = 0; j < NZ; j++) {
                assert(fread(&header,sizeof(int),1,stdin) == 1);
                for (k = 0; k < NX*NY; k++) {
                    assert(fread(&tempfloat,sizeof(float),1,stdin) == 1);
                    ts->dp[l].pos[i] = tempfloat;
                    l++;
                    }
                assert(fread(&trailer,sizeof(int),1,stdin) == 1);
                assert(header == trailer);
                }
	    }
	}
    else if (positionprecision == 1) {
        ts = NULL;
        tsdpp = malloc(sizeof(TIPSY_STRUCTURE));
        assert(tsdpp != NULL);
        tsdpp->th = &th;
        tsdpp->gpdpp = NULL;
        tsdpp->dpdpp = malloc(NX*NY*NZ*sizeof(DARK_PARTICLE));
        assert(tsdpp->dpdpp != NULL);
        tsdpp->spdpp = NULL;
        for (i = 0; i < 3; i++) {
            l = 0;
            for(j = 0; j < NZ; j++) {
                assert(fread(&header,sizeof(int),1,stdin) == 1);
                for (k = 0; k < NX*NY; k++) {
                    assert(fread(&tempfloat,sizeof(float),1,stdin) == 1);
                    tsdpp->dpdpp[l].vel[i] = tempfloat;
                    l++;
                    }
                assert(fread(&trailer,sizeof(int),1,stdin) == 1);
                assert(header == trailer);
                }
            l = 0;
            for(j = 0; j < NZ; j++) {
                assert(fread(&header,sizeof(int),1,stdin) == 1);
                for (k = 0; k < NX*NY; k++) {
                    assert(fread(&tempfloat,sizeof(float),1,stdin) == 1);
                    tsdpp->dpdpp[l].pos[i] = tempfloat;
                    l++;
                    }
                assert(fread(&trailer,sizeof(int),1,stdin) == 1);
                assert(header == trailer);
                }
            }
        }

    /*
    ** Transform coordinates to tipsy unit system and set other properties
    */

    if (particlesoftening == 0) {
	particlesoftening = 1/(NX*20);
	}
    particlemass = (OmB+OmX)/th.ntotal;
    posscalefac = 1/(NX*cell);
    TU_Tipsy = (H_Tipsy*1000)/(h100*100); /* kpc*km^-1*s */
    LU_Tipsy = NX*cell*1000/h100; /* kpc */
    VU_Tipsy = LU_Tipsy/TU_Tipsy; /* km s^-1 */
    velscalefac = 1/(aBegin*VU_Tipsy);
    MU_Tipsy = rhocrit*h100*h100*LU_Tipsy*LU_Tipsy*LU_Tipsy;

    if (positionprecision == 0) {
	for (i = 0; i < th.ndark; i++) {
	    ts->dp[i].mass = particlemass;
	    ts->dp[i].eps = particlesoftening;
	    ts->dp[i].phi = 0;
	    for (j = 0; j < 3; j++) {
		ts->dp[i].pos[j] *= posscalefac;
		ts->dp[i].pos[j] += dr[j];
		ts->dp[i].vel[j] *= velscalefac;
		}
	    }
	}
    else if (positionprecision == 1) {
        for (i = 0; i < th.ndark; i++) {
            tsdpp->dpdpp[i].mass = particlemass;
            tsdpp->dpdpp[i].eps = particlesoftening;
            tsdpp->dpdpp[i].phi = 0;
            for (j = 0; j < 3; j++) {
                tsdpp->dpdpp[i].pos[j] *= posscalefac;
		tsdpp->dpdpp[i].pos[j] += dr[j];
                tsdpp->dpdpp[i].vel[j] *= velscalefac;
                }
            }
	}

    /*
    ** Write out tipsy standard file format
    */

    if (positionprecision == 0) {
	write_tipsy_standard(stdout,ts);
	}
    else if (positionprecision == 1) {
	write_tipsy_standard_dpp(stdout,tsdpp);
	}

    /*
    ** Write out some additional stuff depending on verbose level
    */

    if (verboselevel >= 1) {
	fprintf(stderr,"Paramters from general initial conditions file:\n");
	fprintf(stderr,"name    : %s\n",name);
	fprintf(stderr,"OmB     : %.6e\n",OmB);
	fprintf(stderr,"OmX     : %.6e\n",OmX);
	fprintf(stderr,"OmL     : %.6e\n",OmL);
	fprintf(stderr,"OmN     : %.6e\n",OmN);
	fprintf(stderr,"h100    : %.6e\n",h100);
	fprintf(stderr,"cell    : %.6e chimp\n",cell);
	fprintf(stderr,"enn     : %.6e\n",enn);
	fprintf(stderr,"sigma8  : %.6e\n",sigma8);
	fprintf(stderr,"akpivot : %.6e\n",akpivot);
	fprintf(stderr,"aBegin  : %.6e\n",aBegin);
	fprintf(stderr,"delDC   : %.6e\n",delDC);
	fprintf(stderr,"NX      : %d\n",NX);
	fprintf(stderr,"NY      : %d\n",NY);
	fprintf(stderr,"NZ      : %d\n",NZ);
	fprintf(stderr,"Used values:\n");
	fprintf(stderr,"drx  : %.6e LU\n",dr[0]);
	fprintf(stderr,"dry  : %.6e LU\n",dr[1]);
	fprintf(stderr,"drz  : %.6e LU\n",dr[2]);
	fprintf(stderr,"soft : %.6e LU\n",particlesoftening);
	fprintf(stderr,"Resulting internal tipsy units:\n");
	fprintf(stderr,"LU : %.6e kpc\n",LU_Tipsy);
	fprintf(stderr,"TU : %.6e Gyr\n",TU_Tipsy/VelConvertFac);
	fprintf(stderr,"VU : %.6e km s^-1 = %.6e kpc Gyr^-1\n",VU_Tipsy,VU_Tipsy*VelConvertFac);
	fprintf(stderr,"MU : %.6e Mo\n",MU_Tipsy);
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
    fprintf(stderr,"-spp          : set this flag if input and output files have single precision positions (default)\n");
    fprintf(stderr,"-dpp          : set this flag if input and output files have double precision positions\n");
    fprintf(stderr,"-drx <value>  : shift along x-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-dry <value>  : shift along y-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-drz <value>  : shift along z-axis [LU] (default: -0.5 LU)\n");
    fprintf(stderr,"-soft <value> : softening length of particles [LU] (default: 1/20 mean particle separation => 1/(N^{-3}*20) LU)\n");
    fprintf(stderr,"-v            : more informative output to screen\n");
    fprintf(stderr,"< <name>      : input file in tipsy binary format\n");
    fprintf(stderr,"> <name>      : output file in tipsy standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

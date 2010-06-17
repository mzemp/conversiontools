/* 
** art_nb_2_silo.c
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
#include <silo.h>
#include <iof.h>
#include <art_sfc.h>

void usage(void);

int main(int argc, char **argv) {

    int positionprecision, verboselevel, scaling, massdarkfromdata, Lmaxgaswrite;
    int L;
    int index[3] = {-1,-1,-1};
    int *cellrefined = NULL;
    long int i, j, k;
    long int mothercellindex, childcellindex;
    long int Nparticleread, Nrecordread, Ngasread, Ngaswritten, Ndarkread, Nstarread;
    long int *Icoordinates = NULL;
    double ***coordinates = NULL;
    double r[3];
    double celllength, cellvolume;
    double softfac, LBox;
    COSMOLOGICAL_PARAMETERS cp;
    UNIT_SYSTEM artus, cosmous;
    COORDINATE_TRANSFORMATION art2cosmo_ct;
    ART_DATA ad;
    ART_GAS_PROPERTIES agp;
    ART_STAR_PROPERTIES asp;
    ART_COORDINATES *ac = NULL;
    DBfile *dbfile = NULL;
    char outputname[256];
    int nvar = 8;
    int types[] = {DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR,
		    DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR,
		    DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR};
/*
    char *gas_names[] = {"gas/gas_rx","gas/gas_ry","gas/gas_rz",
			   "gas/gas_vx","gas/gas_vy","gas/gas_vz",
			   "gas/gas_mag_r","gas/gas_mag_v"};
    char *gas_defs[] = {"coord(<gas/gas_pos>)[0]","coord(<gas/gas_pos>)[1]","coord(<gas/gas_pos>)[2]",
			  "dot(<gas/gas_vel>,{1,0,0})","dot(<gas/gas_vel>,{0,1,0})","dot(<gas/gas_vel>,{0,0,1})",
			  "polar_radius(<gas/gas_pos>)","magnitude(<gas/gas_vel>)"};
*/
    char *dark_names[] = {"dark/dark_rx","dark/dark_ry","dark/dark_rz",
			   "dark/dark_vx","dark/dark_vy","dark/dark_vz",
			   "dark/dark_mag_r","dark/dark_mag_v"};
    char *dark_defs[] = {"coord(<dark/dark_pos>)[0]","coord(<dark/dark_pos>)[1]","coord(<dark/dark_pos>)[2]",
			  "dot(<dark/dark_vel>,{1,0,0})","dot(<dark/dark_vel>,{0,1,0})","dot(<dark/dark_vel>,{0,0,1})",
			  "polar_radius(<dark/dark_pos>)","magnitude(<dark/dark_vel>)"};
    char *star_names[] = {"star/star_rx","star/star_ry","star/star_rz",
			   "star/star_vx","star/star_vy","star/star_vz",
			   "star/star_mag_r","star/star_mag_v"};
    char *star_defs[] = {"coord(<star/star_pos>)[0]","coord(<star/star_pos>)[1]","coord(<star/star_pos>)[2]",
			  "dot(<star/star_vel>,{1,0,0})","dot(<star/star_vel>,{0,1,0})","dot(<star/star_vel>,{0,0,1})",
			  "polar_radius(<star/star_pos>)","magnitude(<star/star_vel>)"};
    float **darkposf = NULL, **starposf = NULL;
    double **darkposd = NULL, **starposd = NULL;
    float **darkvel = NULL, **starvel = NULL;
    float *darkmass = NULL, *starmass = NULL;

    /*
    ** Set some default values
    */

    set_default_values_art_data(&ad);
    set_default_values_coordinate_transformation(&art2cosmo_ct);

    artus.LBox = 0;
    artus.Hubble0 = 0;
    artus.rhocrit0 = 0;

    cosmous.LBox = 0;
    cosmous.Hubble0 = 0;
    cosmous.rhocrit0 = 0;

    positionprecision = 0;
    verboselevel = 0;
    scaling = 0;
    massdarkfromdata = 1;
    Lmaxgaswrite = -1;
    softfac = 50;
    Ngaswritten = 0;
    L = 0;

    /*
    ** Read in options
    */

    i = 1;
    while (i < argc) {
        if (strcmp(argv[i],"-spp") == 0) {
            positionprecision = 0;
            i++;
            }
/*
        else if (strcmp(argv[i],"-dpp") == 0) {
            positionprecision = 1;
            i++;
            }
	else if (strcmp(argv[i],"-massdarkfromdata") == 0) {
	    massdarkfromdata = 1;
            i++;
	    }
	else if (strcmp(argv[i],"-massdarknotfromdata") == 0) {
	    massdarkfromdata = 0;
            i++;
	    }
        else if (strcmp(argv[i],"-toplevelsoftdark") == 0) {
            i++;
            if (i >= argc) usage();
            ad.toplevelsoftdark = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-toplevelmassdark") == 0) {
            i++;
            if (i >= argc) usage();
            ad.toplevelmassdark = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-softfac") == 0) {
            i++;
            if (i >= argc) usage();
            softfac = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Lmaxgaswrite") == 0) {
            i++;
            if (i >= argc) usage();
            Lmaxgaswrite = atoi(argv[i]);
            i++;
            }
*/
        else if (strcmp(argv[i],"-LBox") == 0) {
            i++;
            if (i >= argc) usage();
            LBox = atof(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-GRAVITY") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.GRAVITY = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-HYDRO") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.HYDRO = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-ADVECT_SPECIES") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.ADVECT_SPECIES = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-STARFORM") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.STARFORM = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-ENRICH") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.ENRICH = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-ENRICH_SNIa") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.ENRICH_SNIa = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-RADIATIVE_TRANSFER") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.RADIATIVE_TRANSFER = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-ELECTRON_ION_NONEQUILIBRIUM") == 0) {
	    i++;
            if (i >= argc) usage();
	    ad.ELECTRON_ION_NONEQUILIBRIUM = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-headerfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(ad.HeaderFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-coordinatesdatafile") == 0) {
	    ad.darkcontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.CoordinatesDataFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-starpropertiesfile") == 0) {
	    ad.starcontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.StarPropertiesFileName,argv[i]);
            i++;
            }
/*
        else if (strcmp(argv[i],"-gasfile") == 0) {
	    ad.gascontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.GasFileName,argv[i]);
            i++;
            }
*/
	else if (strcmp(argv[i],"-o") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(outputname,argv[i]);
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
    ** Read header file
    */

    prepare_art_data(&ad);

    if (Lmaxgaswrite == -1) Lmaxgaswrite = ad.Lmaxgas;
    assert(Lmaxgaswrite >= 0);
    assert(Lmaxgaswrite <= ad.Lmaxgas);

    cp.OmegaM0 = ad.ah.OmM0;
    cp.OmegaB0 = ad.ah.OmB0;
    cp.OmegaDM0 = cp.OmegaM0 - cp.OmegaB0;
    cp.OmegaL0 = ad.ah.OmL0;
    cp.OmegaK0 = ad.ah.OmK0;
    cp.OmegaR0 = 0;
    cp.h0_100 = ad.ah.h100;

    if(artus.LBox == 0) artus.LBox = ad.ah.Ngrid;
    if(artus.Hubble0 == 0) artus.Hubble0 = 2.0/sqrt(cp.OmegaM0);
    if(artus.rhocrit0 == 0) artus.rhocrit0 = 1/cp.OmegaM0;

    if(cosmous.LBox == 0) cosmous.LBox = LBox;
    if(cosmous.Hubble0 == 0) cosmous.Hubble0 = 100*cp.h0_100*ConversionFactors.km_per_s_2_kpc_per_Gyr/1e3;
    if(cosmous.rhocrit0 == 0) cosmous.rhocrit0 = PhysicalConstants.rho_crit_Cosmology*pow(cp.h0_100,2);

    calculate_units_transformation(artus,cosmous,&art2cosmo_ct);

    /*
    ** Set masses and softenings for dark matter in ART units
    */

    if (ad.toplevelsoftdark == -1) ad.toplevelsoftdark = ad.rootcelllength/softfac;
    if (ad.toplevelmassdark == -1) {
	if (massdarkfromdata == 1) ad.toplevelmassdark = ad.ah.mass[ad.Lmaxdark];
	else ad.toplevelmassdark = cp.OmegaDM0/cp.OmegaM0;
	}
    assert(ad.toplevelsoftdark > 0);
    assert(ad.toplevelmassdark > 0);
    ad.toplevelsoftdark = 0;
    for (i = ad.Lmindark; i <= ad.Lmaxdark; i++) {
	ad.softdark[i] = ad.toplevelsoftdark/pow(ad.refinementstepdark,i);
	if (massdarkfromdata == 1) ad.massdark[i] = ad.ah.mass[ad.Lmaxdark-i];
	else ad.massdark[i] = ad.toplevelmassdark/pow(ad.refinementstepdark,3*i);
	}

    /*
    ** Get silo arrays ready
    */

    if (positionprecision == 0) {
	darkposf = realloc(starposf,3*sizeof(float *));
	assert(darkposf != NULL);
	starposf = realloc(starposf,3*sizeof(float *));
	assert(starposf != NULL);
	}
    else if (positionprecision == 1) {
	darkposd = realloc(darkposd,3*sizeof(double *));
	assert(darkposd != NULL);
	starposd = realloc(starposd,3*sizeof(double *));
	assert(starposd != NULL);
	}
    darkvel = malloc(3*sizeof(float *));
    assert(darkvel != NULL);
    starvel = malloc(3*sizeof(float *));
    assert(starvel != NULL);

    /*
    ** Read and process data
    */

    if (0) {
	/*
	** Gas
	*/
	fprintf(stderr,"Processing gas ... ");
	coordinates = malloc((ad.Lmaxgas+1)*sizeof(double **));
	assert(coordinates != NULL);
	Icoordinates = malloc((ad.Lmaxgas+1)*sizeof(long int));
	assert(Icoordinates != NULL);
	for (i = 0; i < (ad.Lmaxgas+1); i++) {
	    Icoordinates[i] = 0;
	    }
	Ngasread = 0;
	init_sfc(&ad.asfci);
	/*
	** Go through all levels
	*/
	for (i = ad.Lmingas; i <= Lmaxgaswrite; i++) {
	    /*
	    ** Calculate level properties and read level header
	    */
	    celllength = ad.rootcelllength/pow(2,i);
	    cellvolume = celllength*celllength*celllength;
	    read_art_nb_gas_header_level(&ad,i,&cellrefined);
	    /*
	    ** get coordinates array ready
	    */
	    if (i < Lmaxgaswrite) {
		coordinates[i] = malloc(ad.Ncellrefined[i]*sizeof(double *));
		assert(coordinates[i] != NULL);
		for (j = 0; j < ad.Ncellrefined[i]; j++) {
		    coordinates[i][j] = malloc(3*sizeof(double));
		    assert(coordinates[i][j] != NULL);
		    }
		}
	    /*
	    ** Move file positions
	    */
	    move_art_nb_gas_filepositions_level_begin(ad,i);
	    /*
	    ** Go through cells in this level
	    */
	    for (j = 0; j < ad.Ncell[i]; j++) {
		read_art_nb_gas_properties(ad,&agp);
		Ngasread++;
		/*
		** Calculate coordinates
		*/
		if (i == ad.Lmingas) {
		    sfc_coords(ad.asfci,j,index);
		    for (k = 0; k < 3; k++) {
			r[k] = index[k] + 0.5;
			}
		    }
		else {
		    for (k = 0; k < 3; k++) {
			mothercellindex = j/8;
			childcellindex = j%8;
			r[k] = coordinates[i-1][mothercellindex][k] + celllength*art_cell_delta[childcellindex][k];
			}
		    }
		/*
		** Check if cell is refined
		*/
		if ((cellrefined[j] == 0) || (i == Lmaxgaswrite)) {
		    /*
		    ** not refined or maximum level reached => write it out
		    */
		    Ngaswritten++;
		    if (positionprecision == 0) {

			}
		    else if (positionprecision == 1) {

			}
		    }
		else if (i < Lmaxgaswrite) {
		    /*
		    ** refined and lower level than Lmaxgaswrite => add it to corresponding coordinates array
		    */
		    for (k = 0; k < 3; k++) {
			coordinates[i][Icoordinates[i]][k] = r[k];
			}
		    Icoordinates[i]++;
		    }
		}
	    /*
	    ** Move file positions
	    */
	    move_art_nb_gas_filepositions_level_end(ad,i);
	    /*
	    ** Checks and free coordinates of level below
	    */
	    if (i < Lmaxgaswrite) assert(Icoordinates[i] == ad.Ncellrefined[i]);
	    if (i > ad.Lmingas) {
		for (j = 0; j < ad.Ncellrefined[i-1]; j++) {
		    free(coordinates[i-1][j]);
		    }
		free(coordinates[i-1]);
		}
	    }
	/*
	** Some checks and free remaining arrays
	*/
	if (Lmaxgaswrite == ad.Lmaxgas) {
	    assert(ad.Ncellrefined[ad.Lmaxgas] == 0);
	    assert(ad.Ngas == Ngasread);
	    j = 0;
	    k = 0;
	    for (i = ad.Lmingas; i <= ad.Lmaxgas; i++) {
		j += ad.Ncell[i];
		k += ad.Ncellrefined[i];
		}
	    assert(ad.Ngas == j);
	    assert(ad.Ngas == k + Ngaswritten);
	    }
	free(Icoordinates);
	free(cellrefined);
	fprintf(stderr,"Done. Processed in total %ld gas particles whereof %ld written out.\n\n",ad.Ngas,Ngaswritten);
	}
    if (ad.darkcontained || ad.starcontained) {
	/*
	** Dark Matter and Stars
	*/
	fprintf(stderr,"Processing dark matter and stars ... ");
	/*
	** Get silo arrays ready
	*/
	for (j = 0; j < 3; j++) {
	    if (positionprecision == 0) {
		darkposf[j] = realloc(darkposf[j],ad.Ndark*sizeof(float));
		assert(darkposf[j] != NULL);
		starposf[j] = realloc(starposf[j],ad.Nstar*sizeof(float));
		assert(starposf[j] != NULL);
		}
	    else if (positionprecision == 1) {
		darkposd[j] = realloc(darkposd[j],ad.Ndark*sizeof(double));
		assert(darkposd[j] != NULL);
		starposd[j] = realloc(starposd[j],ad.Nstar*sizeof(double));
		assert(starposd[j] != NULL);
		}
	    darkvel[j] = realloc(darkvel[j],ad.Ndark*sizeof(float));
	    assert(darkvel[j] != NULL);
	    starvel[j] = realloc(starvel[j],ad.Nstar*sizeof(float));
	    assert(starvel[j] != NULL);
	    }
	darkmass = realloc(darkmass,ad.Ndark*sizeof(float));
	assert(darkmass != NULL);
	starmass = realloc(starmass,ad.Ndark*sizeof(float));
	assert(starmass != NULL);
	/*
	** Get ART stuff ready
	*/
	ac = malloc(ad.Nparticleperrecord*sizeof(ART_COORDINATES));
	assert(ac != NULL);
	if (ad.starcontained) move_art_nb_star_filepositions_begin(ad);
	Nparticleread = 0;
	Ndarkread = 0;
	Nstarread = 0;
	Nrecordread = 0;
	for (i = 0; i < ad.Nrecord; i++) {
	    read_art_nb_coordinates_record(ad,ac);
	    for (j = 0; j < ad.Nparticleperrecord; j++) {
		if (Nparticleread < ad.Ndark) {
		    /*
		    ** Dark Matter
		    */
		    if (positionprecision == 0) {
			for (k = 0; k < 3; k++) {
			    darkposf[k][Ndarkread] = ac[j].r[k];
			    }
			}
		    else if (positionprecision == 1) {
			for (k = 0; k < 3; k++) {
			    darkposd[k][Ndarkread] = ac[j].r[k];
			    }
			}
		    for (k = 0; k < 3; k++) {
			darkvel[k][Ndarkread] = ac[j].v[k];
			}
		    for (k = ad.Lmaxdark; k >=0; k--) {
			if (ad.ah.num[k] >= Nparticleread) L = ad.Lmaxdark-k;
			}
		    darkmass[Ndarkread] = ad.massdark[L];
		    Nparticleread++;
		    Ndarkread++;
		    }
		else if (Nparticleread < ad.Ndark+ad.Nstar) {
		    /*
		    ** Star
		    */
		    read_art_nb_star_properties(ad,&asp);
		    if (positionprecision == 0) {
			for (k = 0; k < 3; k++) {
			    starposf[k][Nstarread] = ac[j].r[k];
			    }
			}
		    else if (positionprecision == 0) {
			for (k = 0; k < 3; k++) {
			    starposd[k][Nstarread] = ac[j].r[k];
			    }
			}
		    for (k = 0; k < 3; k++) {
			starvel[k][Nstarread] = ac[j].v[k];
			}
		    starmass[Nstarread] = asp.mass;
		    Nparticleread++;
		    Nstarread++;
		    }
		}
	    }
	if (ad.starcontained) move_art_nb_star_filepositions_end(ad);
	free(ac);
	fprintf(stderr,"Done. Processed in total %ld dark matter and %ld star particles.\n\n",ad.Ndark,ad.Nstar);
	}

    fprintf(stderr,"Writing data ... ");

    /*
    ** Create silo file
    */

    dbfile = DBCreate(outputname,DB_CLOBBER,DB_LOCAL,"N-body",DB_PDB);
    assert(dbfile != NULL);

    /*
    ** Write header
    */

    L = 1;
    DBMkDir(dbfile,"header");
    DBSetDir(dbfile,"/header");
    DBWrite(dbfile,"Ndark",&ad.Ndark,&L,1,DB_INT);
    DBWrite(dbfile,"Nstar",&ad.Nstar,&L,1,DB_INT);
    DBWrite(dbfile,"auni",&ad.ah.aunin,&L,1,DB_DOUBLE);
    DBWrite(dbfile,"abox",&ad.ah.abox,&L,1,DB_DOUBLE);
    DBSetDir(dbfile,"/");

    DBMkDir(dbfile,"dark");
    DBSetDir(dbfile,"/dark");
    if (positionprecision == 0) {
	DBPutPointmesh(dbfile,"dark_pos",3,darkposf,ad.Ndark,DB_FLOAT,NULL);
	}
    else if (positionprecision == 1) {
	DBPutPointmesh(dbfile,"dark_pos",3,(float **)darkposd,ad.Ndark,DB_DOUBLE,NULL);
	}
    DBPutPointvar(dbfile,"dark_vel","dark_pos",3,darkvel,ad.Ndark,DB_FLOAT,NULL);
    DBPutPointvar1(dbfile,"dark_mass","dark_pos",darkmass,ad.Ndark,DB_FLOAT,NULL);
    DBPutDefvars(dbfile,"variables",nvar,dark_names,types,dark_defs,NULL);
    DBSetDir(dbfile,"/");

    DBMkDir(dbfile,"star");
    DBSetDir(dbfile,"/star");
    if (positionprecision == 0) {
	DBPutPointmesh(dbfile,"star_pos",3,starposf,ad.Nstar,DB_FLOAT,NULL);
	}
    else if (positionprecision == 1) {
	DBPutPointmesh(dbfile,"star_pos",3,(float **)starposd,ad.Nstar,DB_DOUBLE,NULL);
	}
    DBPutPointvar(dbfile,"star_vel","star_pos",3,starvel,ad.Nstar,DB_FLOAT,NULL);
    DBPutPointvar1(dbfile,"star_mass","star_pos",starmass,ad.Nstar,DB_FLOAT,NULL);
    DBPutDefvars(dbfile,"variables",nvar,star_names,types,star_defs,NULL);
    DBSetDir(dbfile,"/");

    DBClose(dbfile);

    fprintf(stderr,"Done.\n\n");

    /*
    ** Write out some additional stuff depending on verbose level
    */

    if (verboselevel >= 1) {
	fprintf(stderr,"There are %d refinement levels for the dark matter:\n\n",ad.Nleveldark);
	for (L = ad.Lmindark; L <= ad.Lmaxdark; L++) {
	    fprintf(stderr,"L %d Lmax %d Nlevel %d Softening %.6e LU_ART = %.6e kpc Mass %.6e MU_ART = %.6e Mo\n",L,ad.Lmaxdark,ad.Ndarklevel[L],ad.softdark[L],ad.softdark[L]*art2cosmo_ct.L_usf,ad.massdark[L],ad.massdark[L]*art2cosmo_ct.M_usf);
	    }
	fprintf(stderr,"\n");
	fprintf(stderr,"ART general header:\n\n");
        fprintf(stderr,"aunin    : %.6e\n",ad.ah.aunin);
        fprintf(stderr,"auni0    : %.6e\n",ad.ah.auni0);
        fprintf(stderr,"amplt    : %.6e\n",ad.ah.amplt);
        fprintf(stderr,"astep    : %.6e\n",ad.ah.astep);
        fprintf(stderr,"istep    : %d\n",ad.ah.istep);
        fprintf(stderr,"partw    : %.6e\n",ad.ah.partw);
        fprintf(stderr,"tintg    : %.6e\n",ad.ah.tintg);
        fprintf(stderr,"ekin     : %.6e\n",ad.ah.ekin);
        fprintf(stderr,"ekin1    : %.6e\n",ad.ah.ekin1);
        fprintf(stderr,"ekin2    : %.6e\n",ad.ah.ekin2);
        fprintf(stderr,"au0      : %.6e\n",ad.ah.au0);
        fprintf(stderr,"aeu0     : %.6e\n",ad.ah.aeu0);
        fprintf(stderr,"Nrow     : %d\n",ad.ah.Nrow);
        fprintf(stderr,"Ngrid    : %d\n",ad.ah.Ngrid);
        fprintf(stderr,"Nspecies : %d\n",ad.ah.Nspecies);
        fprintf(stderr,"Nseed    : %d\n",ad.ah.Nseed);
        fprintf(stderr,"OmM0     : %.6e\n",ad.ah.OmM0);
        fprintf(stderr,"OmL0     : %.6e\n",ad.ah.OmL0);
        fprintf(stderr,"h100     : %.6e\n",ad.ah.h100);
        fprintf(stderr,"Wp5      : %.6e\n",ad.ah.Wp5);
        fprintf(stderr,"OmK0     : %.6e\n",ad.ah.OmK0);
        fprintf(stderr,"OmB0     : %.6e\n",ad.ah.OmB0);
        fprintf(stderr,"magic1   : %.6e\n",ad.ah.magic1);
        fprintf(stderr,"DelDC    : %.6e\n",ad.ah.DelDC);
        fprintf(stderr,"abox     : %.6e\n",ad.ah.abox);
        fprintf(stderr,"Hbox     : %.6e\n",ad.ah.Hbox);
        fprintf(stderr,"magic2   : %.6e\n",ad.ah.magic2);
        fprintf(stderr,"Banner   : %s\n",ad.Banner);
	for (i = 0; i < 10; i++) {
       	    fprintf(stderr,"mass[%ld] : %.6e num[%ld] : %d\n",i,ad.ah.mass[i],i,ad.ah.num[i]);
	    }
	fprintf(stderr,"\n");
	fprintf(stderr,"ART data properties:\n\n");
        fprintf(stderr,"Nparticleperrecord : %d\n",ad.Nparticleperrecord);
	fprintf(stderr,"Nrecord            : %d\n",ad.Nrecord);
	fprintf(stderr,"Nhydroproperties   : %d\n",ad.Nhydroproperties);
	fprintf(stderr,"Notherproperties   : %d\n",ad.Notherproperties);
	fprintf(stderr,"Nrtchemspecies     : %d\n",ad.Nrtchemspecies);
	fprintf(stderr,"Nchemspecies       : %d\n",ad.Nchemspecies);
	fprintf(stderr,"Nstarproperties    : %d\n",ad.Nstarproperties);
/*
	fprintf(stderr,"Lmingas            : %d\n",ad.Lmingas);
	fprintf(stderr,"Lmaxgas            : %d\n",ad.Lmaxgas);
	fprintf(stderr,"Lmindark           : %d\n",ad.Lmindark);
	fprintf(stderr,"Lmaxdark           : %d\n",ad.Lmaxdark);
	fprintf(stderr,"Toplevelsoftdark   : %.6e LU_ART = %.6e kpc\n",ad.toplevelsoftdark,ad.toplevelsoftdark*art2cosmo_ct.L_usf);
	fprintf(stderr,"Toplevelmassdark   : %.6e MU_ART = %.6e Mo\n",ad.toplevelmassdark,ad.toplevelmassdark*art2cosmo_ct.M_usf);
*/
	fprintf(stderr,"\n");
	fprintf(stderr,"ART preprocessor flags:\n\n");
	fprintf(stderr,"-GRAVITY                     : %s\n",(ad.GRAVITY == 0)?"not set":"set");
	fprintf(stderr,"-HYDRO                       : %s\n",(ad.HYDRO == 0)?"not set":"set");
	fprintf(stderr,"-ADVECT_SPECIES              : %s\n",(ad.ADVECT_SPECIES == 0)?"not set":"set");
	fprintf(stderr,"-STARFORM                    : %s\n",(ad.STARFORM == 0)?"not set":"set");
	fprintf(stderr,"-ENRICH                      : %s\n",(ad.ENRICH == 0)?"not set":"set");
	fprintf(stderr,"-ENRICH_SNIa                 : %s\n",(ad.ENRICH_SNIa == 0)?"not set":"set");
	fprintf(stderr,"-RADIATIVE_TRANSFER          : %s\n",(ad.RADIATIVE_TRANSFER == 0)?"not set":"set");
	fprintf(stderr,"-ELECTRON_ION_NONEQUILIBRIUM : %s\n",(ad.ELECTRON_ION_NONEQUILIBRIUM == 0)?"not set":"set");
	fprintf(stderr,"\n");
	fprintf(stderr,"ART units:\n\n");
	fprintf(stderr,"LU_ART   : %.6e kpc\n",art2cosmo_ct.L_usf);
	fprintf(stderr,"TU_ART   : %.6e Gyr\n",art2cosmo_ct.T_usf);
	fprintf(stderr,"VU_ART   : %.6e kpc Gyr^{-1} = %.6e km s^{-1}\n",art2cosmo_ct.V_usf,art2cosmo_ct.V_usf*ConversionFactors.kpc_per_Gyr_2_km_per_s);
	fprintf(stderr,"MU_ART   : %.6e Mo\n",art2cosmo_ct.M_usf);
	fprintf(stderr,"\n");
        fprintf(stderr,"Cosmology:\n\n");
        fprintf(stderr,"OmegaM0  : %.6e\n",cp.OmegaM0);
        fprintf(stderr,"OmegaDM0 : %.6e\n",cp.OmegaDM0);
        fprintf(stderr,"OmegaB0  : %.6e\n",cp.OmegaB0);
        fprintf(stderr,"OmegaL0  : %.6e\n",cp.OmegaL0);
        fprintf(stderr,"OmegaK0  : %.6e\n",cp.OmegaK0);
        fprintf(stderr,"OmegaR0  : %.6e\n",cp.OmegaR0);
        fprintf(stderr,"h0_100   : %.6e\n",cp.h0_100);
	fprintf(stderr,"\n");
	fprintf(stderr,"Used values:\n\n");
/*
        fprintf(stderr,"softfac                    : %.6e\n",softfac);
        fprintf(stderr,"Lmaxgaswrite               : %d\n",Lmaxgaswrite);
*/
        fprintf(stderr,"LBox                       : %.6e kpc\n",LBox);
/*
        fprintf(stderr,"Position precision         : %s\n",(positionprecision == 0)?"spp":"dpp");
        fprintf(stderr,"Dark matter mass from data : %s\n",(massdarkfromdata == 0)?"no":"yes");
*/
	fprintf(stderr,"\n");
	}
    if (verboselevel >= 0) {
	ad.Ngas = 0;
        fprintf(stderr,"Time: %g Ntotal: %lu Ngas: %lu Ndark: %lu Nstar: %lu\n",
		ad.ah.aunin,ad.Ndark+ad.Nstar,ad.Ngas,ad.Ndark,ad.Nstar);
        }

    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts ART native binary format to silo format.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
/*
    fprintf(stderr,"-spp                                 : set this flag if output file has single precision positions (default)\n");
    fprintf(stderr,"-dpp                                 : set this flag if output file has double precision positions\n");
    fprintf(stderr,"-massdarkfromdata                    : take dark matter particle masses from data (default)\n");
    fprintf(stderr,"-massdarknotfromdata                 : use value of top level dark matter particle mass for calculating dark matter particle masses\n");
    fprintf(stderr,"-toplevelsoftdark <value>            : softening length of top level dark matter particles [LU_ART] (default: rootcelllength/softfac)\n");
    fprintf(stderr,"-toplevelmassdark <value>            : mass of top level dark matter particles [MU_ART] (default: OmegaDM0/OmegaM0)\n");
    fprintf(stderr,"-softfac <value>                     : softening factor (default: 50)\n");
    fprintf(stderr,"-Lmaxgaswrite <value>                : maximum level of gas written out [counting from 0] (default: Lmaxgas in data)\n");
*/
    fprintf(stderr,"-LBox <value>                        : comoving box length [kpc]\n");
    fprintf(stderr,"-GRAVITY <value>                     : 0 = flag not set / 1 = flag set (default: 1)\n");
    fprintf(stderr,"-HYDRO <value>                       : 0 = flag not set / 1 = flag set (default: 1)\n");
    fprintf(stderr,"-ADVECT_SPECIES <value>              : 0 = flag not set / 1 = flag set (default: 1)\n");
    fprintf(stderr,"-STARFORM <value>                    : 0 = flag not set / 1 = flag set (default: 1)\n");
    fprintf(stderr,"-ENRICH <value>                      : 0 = flag not set / 1 = flag set (default: 1)\n");
    fprintf(stderr,"-ENRICH_SNIa <value>                 : 0 = flag not set / 1 = flag set (default: 1)\n");
    fprintf(stderr,"-RADIATIVE_TRANSFER <value>          : 0 = flag not set / 1 = flag set (default: 1)\n");
    fprintf(stderr,"-ELECTRON_ION_NONEQUILIBRIUM <value> : 0 = flag not set / 1 = flag set (default: 0)\n");
    fprintf(stderr,"-headerfile <name>                   : header file in ART native binary format\n");
    fprintf(stderr,"-coordinatesdatafile <name>          : coordinates data file in ART native binary format\n");
    fprintf(stderr,"-starpropertiesfile <name>           : star properties file in ART native binary format\n");
/*
    fprintf(stderr,"-gasfile <name>                      : gas file in ART native binary format\n");
*/
    fprintf(stderr,"-v                                   : more informative output to screen\n");
    fprintf(stderr,"-o <name>                            : output file in silo format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

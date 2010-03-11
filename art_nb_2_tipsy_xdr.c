/* 
** art_nb_2_tipsy_xdr.c
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
#include <art_sfc.h>

void usage(void);

int main(int argc, char **argv) {

    int positionprecision, verboselevel, scaling, massdarkfromdata, Lmaxgaswrite;
    int L;
    int index[3] = {-1,-1,-1};
    int *cellrefined = NULL;
    long int i, j, k;
    long int mothercellindex, childcellindex;
    long int Nparticleread, Nrecordread, Ngasread, Ngaswritten;
    long int *Icoordinates = NULL;
    double ***coordinates = NULL;
    double r[3];
    double celllength, cellvolume;
    double softfac, LBox;
    double chemical_species_number_density;
    COSMOLOGICAL_PARAMETERS cp;
    UNIT_SYSTEM artus, tipsyus, cosmous;
    COORDINATE_TRANSFORMATION art2tipsy_ct, art2cosmo_ct, tipsy2cosmo_ct;
    TIPSY_HEADER th;
    TIPSY_GAS_PARTICLE tgp;
    TIPSY_DARK_PARTICLE tdp;
    TIPSY_STAR_PARTICLE tsp;
    TIPSY_GAS_PARTICLE_DPP tgpdpp;
    TIPSY_DARK_PARTICLE_DPP tdpdpp;
    TIPSY_STAR_PARTICLE_DPP tspdpp;
    ART_DATA ad;
    ART_GAS_PROPERTIES agp;
    ART_STAR_PROPERTIES asp;
    ART_COORDINATES *ac = NULL;
    FILE *file;
    XDR xdrsin, xdrsout;

    /*
    ** Set some default values
    */

    set_default_values_art_data(&ad);
    set_default_values_coordinate_transformation(&art2tipsy_ct);
    set_default_values_coordinate_transformation(&art2cosmo_ct);
    set_default_values_coordinate_transformation(&tipsy2cosmo_ct);

    artus.LBox = 0;
    artus.Hubble0 = 0;
    artus.rhocrit0 = 0;

    tipsyus.LBox = 0;
    tipsyus.Hubble0 = 0;
    tipsyus.rhocrit0 = 0;

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
        else if (strcmp(argv[i],"-dpp") == 0) {
            positionprecision = 1;
            i++;
            }
	else if (strcmp(argv[i],"-scaling") == 0) {
            scaling = 1;
            i++;
	    }
	else if (strcmp(argv[i],"-noscaling") == 0) {
            scaling = 0;
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
        else if (strcmp(argv[i],"-gasfile") == 0) {
	    ad.gascontained = 1;
            i++;
            if (i >= argc) usage();
            strcpy(ad.GasFileName,argv[i]);
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
    ** Read heder file
    */

    prepare_art_data(&ad);

    if (Lmaxgaswrite == -1) Lmaxgaswrite = ad.Lmaxgas;
    assert(Lmaxgaswrite >= 0);

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

    if(tipsyus.LBox == 0) tipsyus.LBox = 1;
    if(tipsyus.Hubble0 == 0) tipsyus.Hubble0 = sqrt(8.0*M_PI/3.0);
    if(tipsyus.rhocrit0 == 0) tipsyus.rhocrit0 = 1;

    if(cosmous.LBox == 0) cosmous.LBox = LBox;
    if(cosmous.Hubble0 == 0) cosmous.Hubble0 = 100*cp.h0_100*ConversionFactors.km_per_s_2_kpc_per_Gyr/1e3;
    if(cosmous.rhocrit0 == 0) cosmous.rhocrit0 = PhysicalConstants.rho_crit_Cosmology*pow(cp.h0_100,2);

    /*
    ** Calculate coordinate transformation
    */

    if (scaling == 1) {
	/* 
	** ART canonical momentum, Tipsy comoving velocity
	*/
	art2tipsy_ct.V_cssf = 1/pow(ad.ah.aunin,2);
	/*
	** Box is shifted by 0.5 LBox in Tipsy units
	*/
	for (i = 0; i < 3; i++ ) {
	    art2tipsy_ct.L_css[i] = -0.5;
	    }
	calculate_units_transformation(artus,tipsyus,&art2tipsy_ct);
	}
    calculate_units_transformation(artus,cosmous,&art2cosmo_ct);
    calculate_units_transformation(tipsyus,cosmous,&tipsy2cosmo_ct);

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
    for (i = ad.Lmindark; i <= ad.Lmaxdark; i++) {
	ad.softdark[i] = ad.toplevelsoftdark/pow(ad.refinementstepdark,i);
	if (massdarkfromdata == 1) ad.massdark[i] = ad.ah.mass[ad.Lmaxdark-i];
	else ad.massdark[i] = ad.toplevelmassdark/pow(ad.refinementstepdark,3*i);
	}

    /*
    ** Get output file ready
    */

    th.time = ad.ah.aunin;
    th.ntotal = 0;
    th.ndim = ad.Ndim;
    th.ngas = 0;
    th.ndark = 0;
    th.nstar = 0;

    file = fopen("temporary_file_you_can_delete","w");
    assert(file != NULL);
    xdrstdio_create(&xdrsout,file,XDR_ENCODE);
    write_tipsy_xdr_header(&xdrsout,&th);

    /*
    ** Read and process data
    */

    if (ad.gascontained) {
	/*
	** Gas
	*/
	fprintf(stderr,"Processing gas ... ");
	coordinates = malloc((ad.Lmaxgas+1)*sizeof(double **));
	assert(coordinates != NULL);
	Icoordinates = malloc((ad.Lmaxgas+1)*sizeof(long));
	assert(Icoordinates != NULL);
	for (i = 0; i <= (ad.Lmaxgas+1); i++) {
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
		    chemical_species_number_density  = agp.HI_number_density  + 2*agp.HII_number_density;
		    chemical_species_number_density += agp.HeI_number_density + 2*agp.HeII_number_density + 3*agp.HeIII_number_density;
		    chemical_species_number_density += agp.H2_number_density;
		    if (positionprecision == 0) {
			for (k = 0; k < 3; k++) {
			    tgp.pos[k] = r[k]*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k];
			    tgp.vel[k] = (agp.momentum[k]/agp.gas_density)*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf;
			    }
			tgp.mass =  cellvolume*agp.gas_density*art2tipsy_ct.M_usf;
			tgp.rho = 0; //agp.gas_density*art2tipsy_ct.M_usf/pow(art2tipsy_ct.L_usf,3);
			tgp.temp = 0; //agp.internal_energy*(agp.gamma-1)/chemical_species_number_density*art2tipsy_ct.M_usf*pow(art2tipsy_ct.V_usf,2); // k_B T
			tgp.hsmooth = 0; //celllength*art2tipsy_ct.L_usf;
			tgp.metals = 0; //agp.metal_density_total;
			tgp.phi = 0; //agp.potential*pow(art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf,2);
			write_tipsy_xdr_gas(&xdrsout,&tgp);
			}
		    else if (positionprecision == 1) {
			for (k = 0; k < 3; k++) {
			    tgpdpp.pos[k] = r[k]*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k];
			    tgpdpp.vel[k] = (agp.momentum[k]/agp.gas_density)*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf;
			    }
			tgpdpp.mass =  cellvolume*agp.gas_density*art2tipsy_ct.M_usf;
			tgpdpp.rho = 0;
			tgpdpp.temp = 0;
			tgpdpp.hsmooth = 0;
			tgpdpp.metals = 0;
			tgpdpp.phi = 0;
			write_tipsy_xdr_gas_dpp(&xdrsout,&tgpdpp);
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
	ac = malloc(ad.Nparticleperrecord*sizeof(ART_COORDINATES));
	assert(ac != NULL);
	if (ad.starcontained) move_art_nb_star_filepositions_begin(ad);
	Nparticleread = 0;
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
			    tdp.pos[k] = (ac[j].r[k]-ad.shift)*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k];
			    tdp.vel[k] = ac[j].v[k]*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf;
			    }
			for (k = ad.Lmaxdark; k >=0; k--) {
			    if (ad.ah.num[k] >= Nparticleread) L = ad.Lmaxdark-k;
			    }
			tdp.mass = ad.massdark[L]*art2tipsy_ct.M_usf;
			tdp.eps = 0; //ad.softdark[L]*art2tipsy_ct.L_usf;
			tdp.phi = 0;
			write_tipsy_xdr_dark(&xdrsout,&tdp);
			}
		    else if (positionprecision == 1) {
			for (k = 0; k < 3; k++) {
			    tdpdpp.pos[k] = (ac[j].r[k]-ad.shift)*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k];
			    tdpdpp.vel[k] = ac[j].v[k]*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf;
			    }
			for (k = ad.Lmaxdark; k >=0; k--) {
			    if (ad.ah.num[k] >= Nparticleread) L = ad.Lmaxdark-k;
			    }
			tdpdpp.mass = ad.massdark[L]*art2tipsy_ct.M_usf;
			tdpdpp.eps = 0;
			tdpdpp.phi = 0;
			write_tipsy_xdr_dark_dpp(&xdrsout,&tdpdpp);
			}
		    Nparticleread++;
		    }
		else if (Nparticleread < ad.Ndark+ad.Nstar) {
		    /*
		    ** Star
		    */
		    read_art_nb_star_properties(ad,&asp);
		    if (positionprecision == 0) {
			for (k = 0; k < 3; k++) {
			    tsp.pos[k] = (ac[j].r[k]-ad.shift)*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k];
			    tsp.vel[k] = ac[j].v[k]*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf;
			    }
			tsp.mass = asp.mass*art2tipsy_ct.M_usf;
			tsp.metals = 0; // asp.metallicity_SNII+asp.metallicity_SNIa;
			tsp.tform = 0; // requires integral
			tsp.eps = 0;
			tsp.phi = 0;
			write_tipsy_xdr_star(&xdrsout,&tsp);
			}
		    else if (positionprecision == 0) {
			for (k = 0; k < 3; k++) {
			    tspdpp.pos[k] = (ac[j].r[k]-ad.shift)*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k];
			    tspdpp.vel[k] = ac[j].v[k]*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf;
			    }
			tspdpp.mass = asp.mass*art2tipsy_ct.M_usf;
			tspdpp.metals = 0;
			tspdpp.tform = 0;
			tspdpp.eps = 0;
			tspdpp.phi = 0;
			write_tipsy_xdr_star_dpp(&xdrsout,&tspdpp);
			}
		    Nparticleread++;
		    }
		}
	    }
	if (ad.starcontained) move_art_nb_star_filepositions_end(ad);
	free(ac);
	fprintf(stderr,"Done. Processed in total %ld dark matter and %ld star particles.\n\n",ad.Ndark,ad.Nstar);
	}

    xdr_destroy(&xdrsout);
    fclose(file);

    fprintf(stderr,"Rewriting data ... ");

    file = fopen("temporary_file_you_can_delete","r");
    xdrstdio_create(&xdrsin,file,XDR_DECODE);
    xdrstdio_create(&xdrsout,stdout,XDR_ENCODE);

    read_tipsy_xdr_header(&xdrsin,&th);
    th.time = ad.ah.aunin;
    th.ntotal = Ngaswritten+ad.Ndark+ad.Nstar;
    th.ndim = ad.Ndim;
    th.ngas = Ngaswritten; 
    th.ndark = ad.Ndark;
    th.nstar = ad.Nstar;
    write_tipsy_xdr_header(&xdrsout,&th);
    if (positionprecision == 0) {
	for (i = 0; i < th.ngas; i++) { 
	    read_tipsy_xdr_gas(&xdrsin,&tgp);
	    write_tipsy_xdr_gas(&xdrsout,&tgp);
	    }
	for (i = 0; i < th.ndark; i++) { 
	    read_tipsy_xdr_dark(&xdrsin,&tdp);
	    write_tipsy_xdr_dark(&xdrsout,&tdp);
	    }
	for (i = 0; i < th.nstar; i++) { 
	    read_tipsy_xdr_star(&xdrsin,&tsp);
	    write_tipsy_xdr_star(&xdrsout,&tsp);
	    }
	}
    xdr_destroy(&xdrsin);
    xdr_destroy(&xdrsout);
    fclose(file);

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
        fprintf(stderr,"Banner   : %s\n",ad.Banner);
	fprintf(stderr,"\n");
	fprintf(stderr,"ART data properties:\n\n");
        fprintf(stderr,"Nparticleperrecord : %d\n",ad.Nparticleperrecord);
	fprintf(stderr,"Nhydroproperties   : %d\n",ad.Nhydroproperties);
	fprintf(stderr,"Notherproperties   : %d\n",ad.Notherproperties);
	fprintf(stderr,"Nrtchemspecies     : %d\n",ad.Nrtchemspecies);
	fprintf(stderr,"Nchemspecies       : %d\n",ad.Nchemspecies);
	fprintf(stderr,"Toplevelsoftdark   : %.6e LU_ART = %.6e kpc\n",ad.toplevelsoftdark,ad.toplevelsoftdark*art2cosmo_ct.L_usf);
	fprintf(stderr,"Toplevelmassdark   : %.6e MU_ART = %.6e Mo\n",ad.toplevelmassdark,ad.toplevelmassdark*art2cosmo_ct.M_usf);
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
	fprintf(stderr,"Tipsy units:\n\n");
	fprintf(stderr,"LU_TIPSY : %.6e kpc\n",tipsy2cosmo_ct.L_usf);
	fprintf(stderr,"TU_TIPSY : %.6e Gyr\n",tipsy2cosmo_ct.T_usf);
	fprintf(stderr,"VU_TIPSY : %.6e kpc Gyr^{-1} = %.6e km s^{-1}\n",tipsy2cosmo_ct.V_usf,tipsy2cosmo_ct.V_usf*ConversionFactors.kpc_per_Gyr_2_km_per_s);
	fprintf(stderr,"MU_TIPSY : %.6e Mo\n",tipsy2cosmo_ct.M_usf);
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
	fprintf(stderr,"Coordinate transformation:\n\n");
        fprintf(stderr,"L_usf    : %.6e\n",art2tipsy_ct.L_usf);
        fprintf(stderr,"T_usf    : %.6e\n",art2tipsy_ct.T_usf);
        fprintf(stderr,"V_usf    : %.6e\n",art2tipsy_ct.V_usf);
        fprintf(stderr,"M_usf    : %.6e\n",art2tipsy_ct.M_usf);
        fprintf(stderr,"L_cssf   : %.6e\n",art2tipsy_ct.L_cssf);
        fprintf(stderr,"V_cssf   : %.6e\n",art2tipsy_ct.V_cssf);
	for (i = 0; i < 3; i++) fprintf(stderr,"L_css[%ld] : %.6e\n",i,art2tipsy_ct.L_css[i]);
	for (i = 0; i < 3; i++) fprintf(stderr,"V_css[%ld] : %.6e\n",i,art2tipsy_ct.V_css[i]);
	fprintf(stderr,"\n");
	fprintf(stderr,"Used values:\n\n");
        fprintf(stderr,"softfac                    : %.6e\n",softfac);
        fprintf(stderr,"Lmaxgaswrite               : %d\n",Lmaxgaswrite);
        fprintf(stderr,"LBox                       : %.6e kpc\n",LBox);
        fprintf(stderr,"Position precision         : %s\n",(positionprecision == 0)?"spp":"dpp");
        fprintf(stderr,"Scaling                    : %s\n",(scaling == 0)?"no":"yes");
        fprintf(stderr,"Dark matter mass from data : %s\n",(massdarkfromdata == 0)?"no":"yes");
	fprintf(stderr,"\n");
	}
    if (verboselevel >= 0) {
        fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
		th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
        }

    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts ART native binary format to tipsy XDR format.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp                                 : set this flag if output file has single precision positions (default)\n");
    fprintf(stderr,"-dpp                                 : set this flag if output file has double precision positions\n");
    fprintf(stderr,"-scaling                             : scaling from ART to tipsy units\n");
    fprintf(stderr,"-noscaling                           : no scaling from ART to tipsy units (default)\n");
    fprintf(stderr,"-massdarkfromdata                    : take dark matter particle masses from data (default)\n");
    fprintf(stderr,"-massdarknotfromdata                 : use value of top level dark matter particle mass for calculating dark matter particle masses\n");
    fprintf(stderr,"-toplevelsoftdark <value>            : softening length of top level dark matter particles [LU_ART] (default: rootcelllength/softfac)\n");
    fprintf(stderr,"-toplevelmassdark <value>            : mass of top level dark matter particles [MU_ART] (default: OmegaDM0/OmegaM0)\n");
    fprintf(stderr,"-softfac <value>                     : softening factor (default: 50)\n");
    fprintf(stderr,"-Lmaxgaswrite <value>                : maximum level of gas written out [counting from 0] (default: Lmaxgas in data)\n");
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
    fprintf(stderr,"-gasfile <name>                      : gas file in ART native binary format\n");
    fprintf(stderr,"-v                                   : more informative output to screen\n");
    fprintf(stderr,"> <name>                             : output file in tipsy XDR format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

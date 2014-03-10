/* 
** art_nb_2_ascii.c
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
int check_selection(double *, double *, double *);

int main(int argc, char **argv) {

	int verboselevel, Lmaxgaswrite;
	/* int L; */
	int selected;
	int index[3] = {-1,-1,-1};
	int *cellrefined = NULL;
	long int i, j, k;
	long int mothercellindex, childcellindex;
	/* long int Nparticleread, Nrecordread; */
	long int Ngasread, Ngaswritten, Ngasselected;
	/* long int Ndarkselected, Nstarselected; */
	long int *Icoordinates = NULL;
	double ***coordinates = NULL;
	double r[3];
	double celllength, cellvolume;
	double number_density, temperature;
	double LBox;
	double rsel[3], dsel[3], bsel[6], bsim[6];
	double fH2, fH2sel, rhosel;
	COSMOLOGICAL_PARAMETERS cp;
	UNIT_SYSTEM artus, cosmous;
	COORDINATE_TRANSFORMATION art2cosmo_ct;
	ART_DATA ad;
	ART_GAS_PROPERTIES agp;
/*	   ART_STAR_PROPERTIES asp; */
/*	   ART_COORDINATES *ac = NULL; */

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

	LBox = 0;
	verboselevel = 0;
	Lmaxgaswrite = -1;
	fH2sel = 0;
	rhosel = 0;
	Ngaswritten = 0;
	Ngasselected = 0;
	/* Ndarkselected = 0; */
	/* Nstarselected = 0; */
	/* L = 0; */
	for (i = 0; i < 3; i++) {
		rsel[i] = 0;
		dsel[i] = 0;
		bsel[i] = 0;
		bsel[i+3] = 0;
		}

	/*
	** Read in options
	*/

	i = 1;
	while (i < argc) {
		if (strcmp(argv[i],"-pfm") == 0) {
			i++;
			if (i >= argc) usage();
			ad.particle_file_mode = atoi(argv[i]);
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
		else if (strcmp(argv[i],"-rxsel") == 0) {
			i++;
			if (i >= argc) usage();
			rsel[0] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rysel") == 0) {
			i++;
			if (i >= argc) usage();
			rsel[1] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rzsel") == 0) {
			i++;
			if (i >= argc) usage();
			rsel[2] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-dxsel") == 0) {
			i++;
			if (i >= argc) usage();
			dsel[0] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-dysel") == 0) {
			i++;
			if (i >= argc) usage();
			dsel[1] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-dzsel") == 0) {
			i++;
			if (i >= argc) usage();
			dsel[2] = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-fH2sel") == 0) {
			i++;
			if (i >= argc) usage();
			fH2sel = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rhosel") == 0) {
			i++;
			if (i >= argc) usage();
			rhosel = atof(argv[i]);
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
/*		   else if (strcmp(argv[i],"-coordinatesdatafile") == 0) { */
/*		ad.darkcontained = 1; */
/*			   i++; */
/*			   if (i >= argc) usage(); */
/*			   strcpy(ad.CoordinatesDataFileName,argv[i]); */
/*			   i++; */
/*			   } */
/*		   else if (strcmp(argv[i],"-starpropertiesfile") == 0) { */
/*		ad.starcontained = 1; */
/*			   i++; */
/*			   if (i >= argc) usage(); */
/*			   strcpy(ad.StarPropertiesFileName,argv[i]); */
/*			   i++; */
/*			   } */
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
			fprintf(stderr,"Input argument error: %s\n",argv[i]);
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
	cp.OmegaD0 = cp.OmegaM0 - cp.OmegaB0;
	cp.OmegaL0 = ad.ah.OmL0;
	cp.OmegaK0 = ad.ah.OmK0;
	cp.OmegaR0 = 0;
	cp.h0_100 = ad.ah.h100;

	if(artus.LBox == 0) artus.LBox = ad.ah.Ngrid;
	if(artus.Hubble0 == 0) artus.Hubble0 = 2.0/sqrt(cp.OmegaM0);
	if(artus.rhocrit0 == 0) artus.rhocrit0 = 1/cp.OmegaM0;

	if(cosmous.LBox == 0) cosmous.LBox = LBox;
	if(cosmous.Hubble0 == 0) cosmous.Hubble0 = 100*cp.h0_100*ConversionFactors.km_per_s_2_kpc_per_Gyr/1e3;
	if(cosmous.rhocrit0 == 0) cosmous.rhocrit0 = PhysicalConstants.rho_crit_cosmology*pow(cp.h0_100,2);

	/*
	** Calculate coordinate transformation
	*/

	calculate_units_transformation(artus,cosmous,&art2cosmo_ct);

	/*
	** Calculate selection boundaries
	*/

	for (i = 0; i < 3; i++) {
		rsel[i] = put_in_box(rsel[i],0,artus.LBox);
		if (dsel[i] > 0) {
			bsim[i] = 0;
			bsim[i+3] = artus.LBox;
			bsel[i] = rsel[i] - 0.5*dsel[i];
			bsel[i+3] = rsel[i] + 0.5*dsel[i];
			}
		}

	/*
	** Get output file ready
	*/

	fprintf(stdout,"#rx/1 ry/2 rz/3 vx/4 vy/5 vz/6 l/7 V/8");
	fprintf(stdout," rho_gas/9 rho_HI/10 rho_HII/11 rho_H2/12 rho_HeI/13 rho_HeII/14 rho_HeIII/15 metal_density_SNII/16 metal_density_SNIa/17");
	fprintf(stdout," energy/18 pressure/19 gamma/20 internal_energy/21 temperature/22");
	fprintf(stdout,"\n");

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
					selected = check_selection(r,bsim,bsel);
					if (selected) {
						/*
						** H2 fraction
						*/
						fH2 = agp.H2_density/(agp.HI_density+agp.HII_density+agp.H2_density); /* mass densities */
						/*
						** Number density 1/LU^3
						*/
						number_density = agp.HI_density*PhysicalConstants.Mo*art2cosmo_ct.M_usf/PhysicalConstants.m_proton;
						number_density += 2*agp.HII_density*PhysicalConstants.Mo*art2cosmo_ct.M_usf/PhysicalConstants.m_proton;
						number_density += agp.HeI_density*PhysicalConstants.Mo*art2cosmo_ct.M_usf/(4*PhysicalConstants.m_proton);
						number_density += 2*agp.HeII_density*PhysicalConstants.Mo*art2cosmo_ct.M_usf/(4*PhysicalConstants.m_proton);
						number_density += 3*agp.HeIII_density*PhysicalConstants.Mo*art2cosmo_ct.M_usf/(4*PhysicalConstants.m_proton);
						number_density += agp.H2_density*PhysicalConstants.Mo*art2cosmo_ct.M_usf/(2*PhysicalConstants.m_proton);
						/*
						** Internal energy => internal energy per unit volume
						*/
						temperature = agp.internal_energy*(agp.gamma-1)/number_density;
						/*
						** Convert to Joules
						*/
						temperature *= art2cosmo_ct.M_usf*pow(art2cosmo_ct.V_usf,2);
						temperature *= PhysicalConstants.Mo*pow(ConversionFactors.kpc_per_Gyr_2_km_per_s*1000,2);
						/*
						** Gas temperature in Kelvin
						*/
						temperature /= PhysicalConstants.k_Boltzmann;
						/*
						** Check selection
						*/
						if (agp.gas_density >= rhosel || fH2 >= fH2sel) {
							Ngasselected++;
							fprintf(stdout,"%5.4e %5.4e %5.4e",r[0],r[1],r[2]);
							fprintf(stdout," %5.4e %5.4e %5.4e",agp.momentum[0]/agp.gas_density,agp.momentum[1]/agp.gas_density,agp.momentum[2]/agp.gas_density);
							fprintf(stdout," %ld %5.4e %5.4e",i,cellvolume,agp.gas_density);
							fprintf(stdout," %5.4e %5.4e %5.4e",agp.HI_density,agp.HII_density,agp.H2_density);
							fprintf(stdout," %5.4e %5.4e %5.4e",agp.HeI_density,agp.HeII_density,agp.HeIII_density);
							fprintf(stdout," %5.4e %5.4e",agp.metal_density_SNII,agp.metal_density_SNIa);
							fprintf(stdout," %5.4e %5.4e %5.4e %5.4e %5.4e",agp.gas_energy,agp.pressure,agp.gamma,agp.internal_energy,temperature);
							fprintf(stdout,"\n");
							}
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
		fprintf(stderr,"Done. Processed in total %ld gas particles whereof %ld written out.\n\n",ad.Ngas,Ngasselected);
		}
/*	   if ((ad.darkcontained || ad.starcontained) && (writedark || writestar)) { */
/*	/\* */
/*	** Dark Matter and Stars */
/*	*\/ */
/*	fprintf(stderr,"Processing dark matter and stars ... "); */
/*	ac = malloc(ad.Nparticleperrecord*sizeof(ART_COORDINATES)); */
/*	assert(ac != NULL); */
/*	if (ad.starcontained) move_art_nb_star_filepositions_begin(ad); */
/*	Nparticleread = 0; */
/*	Nrecordread = 0; */
/*	for (i = 0; i < ad.Nrecord; i++) { */
/*		read_art_nb_coordinates_record(ad,ac); */
/*		for (j = 0; j < ad.Nparticleperrecord; j++) { */
/*		Nparticleread++; */
/*		for (k = 0; k < 3; k++) r[k] = ac[j].r[k]-ad.shift; */
/*		selected = check_selection(r,bsim,bsel); */
/*		if (Nparticleread <= ad.Ndark) { */
/*			/\* */
/*			** Dark Matter */
/*			*\/ */
/*					   if (writedark && selected) { */
/*			Ndarkselected++; */
/*			if (positionprecision == 0) { */
/*				for (k = 0; k < 3; k++) { */
/*				tdp.pos[k] = r[k]*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k]; */
/*				tdp.vel[k] = ac[j].v[k]*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf; */
/*				} */
/*				for (k = ad.Lmaxdark; k >=0; k--) { */
/*				if (ad.ah.num[k] > Nparticleread) L = ad.Lmaxdark-k; */
/*				} */
/*				tdp.mass = ad.massdark[L]*art2tipsy_ct.M_usf; */
/*				tdp.eps = 0; /\* ad.softdark[L]*art2tipsy_ct.L_usf; *\/ */
/*				tdp.phi = 0; */
/*				write_tipsy_xdr_dark(&xdrsout,&tdp); */
/*				} */
/*			else if (positionprecision == 1) { */
/*				for (k = 0; k < 3; k++) { */
/*				tdpdpp.pos[k] = r[k]*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k]; */
/*				tdpdpp.vel[k] = ac[j].v[k]*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf; */
/*				} */
/*				for (k = ad.Lmaxdark; k >=0; k--) { */
/*				if (ad.ah.num[k] >= Nparticleread) L = ad.Lmaxdark-k; */
/*				} */
/*				tdpdpp.mass = ad.massdark[L]*art2tipsy_ct.M_usf; */
/*				tdpdpp.eps = 0; */
/*				tdpdpp.phi = 0; */
/*				write_tipsy_xdr_dark_dpp(&xdrsout,&tdpdpp); */
/*				} */
/*			} */
/*			} */
/*		else if (Nparticleread <= ad.Ndark+ad.Nstar) { */
/*			/\* */
/*			** Star */
/*			*\/ */
/*			read_art_nb_star_properties(ad,&asp); */
/*			if (writestar && selected) { */
/*			Nstarselected++; */
/*			if (positionprecision == 0) { */
/*				for (k = 0; k < 3; k++) { */
/*				tsp.pos[k] = r[k]*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k]; */
/*				tsp.vel[k] = ac[j].v[k]*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf; */
/*				} */
/*				tsp.mass = asp.mass*art2tipsy_ct.M_usf; */
/*				tsp.metals = 0; /\* asp.metallicity_SNII+asp.metallicity_SNIa; *\/ */
/*				tsp.tform = 0; /\* requires integral *\/ */
/*				tsp.eps = 0; */
/*				tsp.phi = 0; */
/*				write_tipsy_xdr_star(&xdrsout,&tsp); */
/*				} */
/*			else if (positionprecision == 0) { */
/*				for (k = 0; k < 3; k++) { */
/*				tspdpp.pos[k] = r[k]*art2tipsy_ct.L_usf + art2tipsy_ct.L_css[k]; */
/*				tspdpp.vel[k] = ac[j].v[k]*art2tipsy_ct.V_usf*art2tipsy_ct.V_cssf; */
/*				} */
/*				tspdpp.mass = asp.mass*art2tipsy_ct.M_usf; */
/*				tspdpp.metals = 0; */
/*				tspdpp.tform = 0; */
/*				tspdpp.eps = 0; */
/*				tspdpp.phi = 0; */
/*				write_tipsy_xdr_star_dpp(&xdrsout,&tspdpp); */
/*				} */
/*			} */
/*			} */
/*		} */
/*		} */
/*	if (ad.starcontained) move_art_nb_star_filepositions_end(ad); */
/*	free(ac); */
/*	fprintf(stderr,"Done. Processed in total %ld dark matter and %ld star particles whereof %ld and %ld written out.\n\n",ad.Ndark,ad.Nstar,Ndarkselected,Nstarselected); */
/*	} */

	/*
	** Write out some additional stuff depending on verbose level
	*/

	if (verboselevel >= 1) {
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
		fprintf(stderr,"Particle File Mode : %d\n",ad.particle_file_mode);
		fprintf(stderr,"Nparticleperrecord : %d\n",ad.Nparticleperrecord);
		fprintf(stderr,"Nrecord            : %d\n",ad.Nrecord);
		fprintf(stderr,"Nhydroproperties   : %d\n",ad.Nhydroproperties);
		fprintf(stderr,"Notherproperties   : %d\n",ad.Notherproperties);
		fprintf(stderr,"Nrtchemspecies     : %d\n",ad.Nrtchemspecies);
		fprintf(stderr,"Nchemspecies       : %d\n",ad.Nchemspecies);
		fprintf(stderr,"Nstarproperties    : %d\n",ad.Nstarproperties);
		fprintf(stderr,"Lmingas            : %d\n",ad.Lmingas);
		fprintf(stderr,"Lmaxgas            : %d\n",ad.Lmaxgas);
		fprintf(stderr,"Lmindark           : %d\n",ad.Lmindark);
		fprintf(stderr,"Lmaxdark           : %d\n",ad.Lmaxdark);
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
		fprintf(stderr,"LU : %.6e kpc\n",art2cosmo_ct.L_usf);
		fprintf(stderr,"TU : %.6e Gyr\n",art2cosmo_ct.T_usf);
		fprintf(stderr,"VU : %.6e kpc Gyr^{-1} = %.6e km s^{-1}\n",art2cosmo_ct.V_usf,art2cosmo_ct.V_usf*ConversionFactors.kpc_per_Gyr_2_km_per_s);
		fprintf(stderr,"MU : %.6e Mo\n",art2cosmo_ct.M_usf);
		fprintf(stderr,"\n");
		fprintf(stderr,"Cosmology:\n\n");
		fprintf(stderr,"OmegaM0 : %.6e\n",cp.OmegaM0);
		fprintf(stderr,"OmegaD0 : %.6e\n",cp.OmegaD0);
		fprintf(stderr,"OmegaB0 : %.6e\n",cp.OmegaB0);
		fprintf(stderr,"OmegaL0 : %.6e\n",cp.OmegaL0);
		fprintf(stderr,"OmegaK0 : %.6e\n",cp.OmegaK0);
		fprintf(stderr,"OmegaR0 : %.6e\n",cp.OmegaR0);
		fprintf(stderr,"h0_100  : %.6e\n",cp.h0_100);
		fprintf(stderr,"\n");
		fprintf(stderr,"Selected volume:\n\n");
		if (dsel[0] == 0 && dsel[1] == 0 && dsel[2] == 0) {
			fprintf(stderr,"No selection chosen.\n");
			fprintf(stderr,"\n");
			}
		else {
			fprintf(stderr,"Centre: (%.6e, %.6e, %.6e) LU\n",rsel[0],rsel[1],rsel[2]);
			fprintf(stderr,"Width: (%.6e, %.6e, %.6e) LU\n",dsel[0],dsel[1],dsel[2]);
			fprintf(stderr,"Box: [%.6e ... %.6e] x [%.6e ... %.6e] x [%.6e ... %.6e] LU\n",bsel[0],bsel[3],bsel[1],bsel[4],bsel[2],bsel[5]);
			fprintf(stderr,"\n");
			}
		fprintf(stderr,"Used values:\n\n");
		fprintf(stderr,"Lmaxgaswrite : %d\n",Lmaxgaswrite);
		fprintf(stderr,"LBox         : %.6e kpc\n",LBox);
		fprintf(stderr,"fH2sel       : %.6e\n",fH2sel);
		fprintf(stderr,"rhosel       : %.6e MU LU^{-3}\n",rhosel);
		fprintf(stderr,"\n");
		}

	exit(0);
	}

void usage(void) {

	fprintf(stderr,"\n");
	fprintf(stderr,"Program converts ART native binary format to ascii text.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Please specify the following parameters:\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"-pfm <value>                         : particle file mode of ART file (default: 0)\n");
	fprintf(stderr,"-Lmaxgaswrite <value>                : maximum level of gas written out [counting from 0] (default: Lmaxgas in data)\n");
	fprintf(stderr,"-rxsel <value>                       : x-coordinate of centre of selection box [LU] (default: 0)\n");
	fprintf(stderr,"-rysel <value>                       : y-coordinate of centre of selection box [LU] (default: 0)\n");
	fprintf(stderr,"-rzsel <value>                       : z-coordinate of centre of selection box [LU] (default: 0)\n");
	fprintf(stderr,"-dxsel <value>                       : x-width of selection box [LU] (default: 0, i.e. no selection)\n");
	fprintf(stderr,"-dysel <value>                       : y-width of selection box [LU] (default: 0, i.e. no selection)\n");
	fprintf(stderr,"-dzsel <value>                       : z-width of selection box [LU] (default: 0, i.e. no selection)\n");
	fprintf(stderr,"-LBox <value>                        : comoving box length [kpc]\n");
	fprintf(stderr,"-fH2sel <value>                      : molecular hydrogen fraction for selection \n");
	fprintf(stderr,"-rhosel <value>                      : comoving density for selection [MU LU^{-3}]\n");
	fprintf(stderr,"-GRAVITY <value>                     : 0 = flag not set / 1 = flag set (default: 1)\n");
	fprintf(stderr,"-HYDRO <value>                       : 0 = flag not set / 1 = flag set (default: 1)\n");
	fprintf(stderr,"-ADVECT_SPECIES <value>              : 0 = flag not set / 1 = flag set (default: 1)\n");
	fprintf(stderr,"-STARFORM <value>                    : 0 = flag not set / 1 = flag set (default: 1)\n");
	fprintf(stderr,"-ENRICH <value>                      : 0 = flag not set / 1 = flag set (default: 1)\n");
	fprintf(stderr,"-ENRICH_SNIa <value>                 : 0 = flag not set / 1 = flag set (default: 1)\n");
	fprintf(stderr,"-RADIATIVE_TRANSFER <value>          : 0 = flag not set / 1 = flag set (default: 1)\n");
	fprintf(stderr,"-ELECTRON_ION_NONEQUILIBRIUM <value> : 0 = flag not set / 1 = flag set (default: 0)\n");
	fprintf(stderr,"-headerfile <name>                   : header file in ART native binary format\n");
/*	   fprintf(stderr,"-coordinatesdatafile <name>          : coordinates data file in ART native binary format\n"); */
/*	   fprintf(stderr,"-starpropertiesfile <name>           : star properties file in ART native binary format\n"); */
	fprintf(stderr,"-gasfile <name>                      : gas file in ART native binary format\n");
	fprintf(stderr,"-v                                   : more informative output to screen\n");
	fprintf(stderr,"> <name>                             : output file in ascii format\n");
	fprintf(stderr,"\n");
	exit(1);
	}

int check_selection(double r[3], double boxsim[6], double boxsel[6]) {

	int i;
	int selected = 1;
	double dsel, LBox;

	for (i = 0; i < 3; i++) {
		dsel = boxsel[i+3]-boxsel[i];
		LBox = boxsim[i+3]-boxsim[i];
		if (dsel > 0) {
			if (boxsim[i] <= boxsel[i] && boxsim[i+3] >= boxsel[i+3]) {
				/*
				** Selection box is contained within simulation box
				*/
				if (r[i] < boxsel[i] || r[i] > boxsel[i+3]) selected = 0;
				}
			else if (boxsel[i] < boxsim[i]) {
				/*
				** Selection box exceeds lower simulation box boundary
				*/
				if (r[i] < boxsel[i]+LBox && r[i] > boxsel[i+3]) selected = 0;
				if (r[i] >= boxsel[i]+LBox) r[i] -= LBox;
				}
			else if (boxsel[i+3] > boxsim[i+3]) {
				/*
				** Selection box exceeds upper simulation box boundary
				*/
				if (r[i] < boxsel[i] && r[i] > boxsel[i+3]-LBox) selected = 0;
				if (r[i] <= boxsel[i+3]-LBox) r[i] += LBox;
				}
			}
		}
	return selected;
	}

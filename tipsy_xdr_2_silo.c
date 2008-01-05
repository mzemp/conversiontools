/* 
** ts2silo.c
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <silo.h>
#include <iof.h>

void usage(void);

int main(int argc, char **argv) {

    int i, j;
    int denfile = 0, psdfile = 0;
    int nvar = 8;
    int types[] = {DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR,
		    DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR,
		    DB_VARTYPE_SCALAR,DB_VARTYPE_SCALAR};
    const char *gas_names[] = {"gas/gas_rx","gas/gas_ry","gas/gas_rz",
			   "gas/gas_vx","gas/gas_vy","gas/gas_vz",
			   "gas/gas_mag_r","gas/gas_mag_v"};
    const char *gas_defs[] = {"coord(<gas/gas_pos>)[0]","coord(<gas/gas_pos>)[1]","coord(<gas/gas_pos>)[2]",
			  "dot(<gas/gas_vel>,{1,0,0})","dot(<gas/gas_vel>,{0,1,0})","dot(<gas/gas_vel>,{0,0,1})",
			  "polar_radius(<gas/gas_pos>)","magnitude(<gas/gas_vel>)"};
    const char *dark_names[] = {"dark/dark_rx","dark/dark_ry","dark/dark_rz",
			   "dark/dark_vx","dark/dark_vy","dark/dark_vz",
			   "dark/dark_mag_r","dark/dark_mag_v"};
    const char *dark_defs[] = {"coord(<dark/dark_pos>)[0]","coord(<dark/dark_pos>)[1]","coord(<dark/dark_pos>)[2]",
			  "dot(<dark/dark_vel>,{1,0,0})","dot(<dark/dark_vel>,{0,1,0})","dot(<dark/dark_vel>,{0,0,1})",
			  "polar_radius(<dark/dark_pos>)","magnitude(<dark/dark_vel>)"};
    const char *star_names[] = {"star/star_rx","star/star_ry","star/star_rz",
			   "star/star_vx","star/star_vy","star/star_vz",
			   "star/star_mag_r","star/star_mag_v"};
    const char *star_defs[] = {"coord(<star/star_pos>)[0]","coord(<star/star_pos>)[1]","coord(<star/star_pos>)[2]",
			  "dot(<star/star_vel>,{1,0,0})","dot(<star/star_vel>,{0,1,0})","dot(<star/star_vel>,{0,0,1})",
			  "polar_radius(<star/star_pos>)","magnitude(<star/star_vel>)"};
    char outputname[100], denfilename[100], psdfilename[100];
    float **pos;
    float **vel;
    float *mass;
    float *eps;
    float *phi;
    float *rho;
    float *temp;
    float *metals;
    float *tform;
    TIPSY_HEADER th;
    GAS_PARTICLE gp;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
    DBfile *dbfile = NULL;
    FILE *file;
    XDR xdrs;

    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-o") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            strcpy(outputname,argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-den") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
	    denfile = 1;
            strcpy(denfilename,argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-psd") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
	    psdfile = 1;
            strcpy(psdfilename,argv[i]);
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
    ** Create silo file
    */
    dbfile = DBCreate(outputname,DB_CLOBBER,DB_LOCAL,"N-body",DB_PDB);
    assert(dbfile != NULL);
    /*
    ** Read in tipsy file, write out header stuff and get arrays ready
    */
    xdrstdio_create(&xdrs,stdin,XDR_DECODE);
    read_tipsy_standard_header(&xdrs,&th);
    i = 1;
    DBMkDir(dbfile,"header");
    DBSetDir(dbfile,"/header");
    DBWrite(dbfile,"ntotal",&th.ntotal,&i,1,DB_INT);
    DBWrite(dbfile,"ngas",&th.ngas,&i,1,DB_INT);
    DBWrite(dbfile,"ndark",&th.ndark,&i,1,DB_INT);
    DBWrite(dbfile,"nstar",&th.nstar,&i,1,DB_INT);
    DBWrite(dbfile,"time",&th.time,&i,1,DB_DOUBLE);
    DBWrite(dbfile,"ndim",&th.ndim,&i,1,DB_INT);
    DBSetDir(dbfile,"/");
    pos = malloc(th.ndim*sizeof(float *));
    assert(pos != NULL);
    vel = malloc(th.ndim*sizeof(float *));
    assert(vel != NULL);
    for (j = 0; j < th.ndim; j++) {
	pos[j] = NULL;
	vel[j] = NULL;
	}
    mass = NULL;
    eps = NULL;
    phi = NULL;
    rho = NULL;
    temp = NULL;
    metals = NULL;
    tform = NULL;
    /*
    ** Process gas particles
    */
    if (th.ngas > 0) {
	for (j = 0; j < th.ndim; j++) {
	    pos[j] = realloc(pos[j],th.ngas*sizeof(float));
	    assert(pos[j] != NULL);
	    vel[j] = realloc(vel[j],th.ngas*sizeof(float));
	    assert(vel[j] != NULL);
	    }
	mass = realloc(mass,th.ngas*sizeof(float));
	assert(mass != NULL);
	eps = realloc(eps,th.ngas*sizeof(float));
	assert(eps != NULL);
	phi = realloc(phi,th.ngas*sizeof(float));
	assert(phi != NULL);
	rho = realloc(rho,th.ngas*sizeof(float));
	assert(rho != NULL);
	temp = realloc(temp,th.ngas*sizeof(float));
	assert(temp != NULL);
	metals = realloc(metals,th.ngas*sizeof(float));
	assert(metals != NULL);
	for (i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas(&xdrs,&gp);
	    for (j = 0; j < th.ndim; j++) {
		pos[j][i] = gp.pos[j];
		vel[j][i] = gp.vel[j];
		}
	    mass[i] = gp.mass;
	    eps[i] = gp.hsmooth;
	    phi[i] = gp.phi;
	    rho[i] = gp.metals;
	    temp[i] = gp.temp;
	    metals[i] = gp.metals;
	    }
	DBMkDir(dbfile,"gas");
	DBSetDir(dbfile,"/gas");
	DBPutPointmesh(dbfile,"gas_pos",th.ndim,pos,th.ngas,DB_FLOAT,NULL);
	DBPutPointvar(dbfile,"gas_vel","gas_pos",th.ndim,vel,th.ngas,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"gas_mass","gas_pos",mass,th.ngas,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"gas_hsmooth","gas_pos",eps,th.ngas,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"gas_phi","gas_pos",phi,th.ngas,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"gas_metals","gas_pos",phi,th.ngas,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"gas_rho","gas_pos",phi,th.ngas,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"gas_temp","gas_pos",phi,th.ngas,DB_FLOAT,NULL);
	DBPutDefvars(dbfile,"variables",nvar,gas_names,types,gas_defs,NULL);
	DBSetDir(dbfile,"/");
	}
    /*
    ** Process dark matter particles
    */
    if (th.ndark > 0) {
	if (th.ngas > 0) {
	    metals = realloc(metals,0*sizeof(float));
	    rho = realloc(rho,0*sizeof(float));
	    temp = realloc(rho,0*sizeof(float));
	    }
	for (j = 0; j < th.ndim; j++) {
	    pos[j] = realloc(pos[j],th.ndark*sizeof(float));
	    assert(pos[j] != NULL);
	    vel[j] = realloc(vel[j],th.ndark*sizeof(float));
	    assert(vel[j] != NULL);
	    }
	mass = realloc(mass,th.ndark*sizeof(float));
	assert(mass != NULL);
	eps = realloc(eps,th.ndark*sizeof(float));
	assert(eps != NULL);
	phi = realloc(phi,th.ndark*sizeof(float));
	assert(phi != NULL);
	for (i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark(&xdrs,&dp);
	    for (j = 0; j < th.ndim; j++) {
		pos[j][i] = dp.pos[j];
		vel[j][i] = dp.vel[j];
		}
	    mass[i] = dp.mass;
	    eps[i] = dp.eps;
	    phi[i] = dp.phi;
	    }
	DBMkDir(dbfile,"dark");
	DBSetDir(dbfile,"/dark");
	DBPutPointmesh(dbfile,"dark_pos",th.ndim,pos,th.ndark,DB_FLOAT,NULL);
	DBPutPointvar(dbfile,"dark_vel","dark_pos",th.ndim,vel,th.ndark,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"dark_mass","dark_pos",mass,th.ndark,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"dark_eps","dark_pos",eps,th.ndark,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"dark_phi","dark_pos",phi,th.ndark,DB_FLOAT,NULL);
	DBPutDefvars(dbfile,"variables",nvar,dark_names,types,dark_defs,NULL);
	DBSetDir(dbfile,"/");
	}
    /*
    ** Process star particles
    */
    if (th.nstar > 0) {
	for (j = 0; j < th.ndim; j++) {
	    pos[j] = realloc(pos[j],th.nstar*sizeof(float));
	    assert(pos[j] != NULL);
	    vel[j] = realloc(vel[j],th.nstar*sizeof(float));
	    assert(vel[j] != NULL);
	    }
	mass = realloc(mass,th.nstar*sizeof(float));
	assert(mass != NULL);
	eps = realloc(eps,th.nstar*sizeof(float));
	assert(eps != NULL);
	phi = realloc(phi,th.nstar*sizeof(float));
	assert(phi != NULL);
	metals = realloc(metals,th.nstar*sizeof(float));
	assert(metals != NULL);
	tform = realloc(tform,th.nstar*sizeof(float));
	assert(tform != NULL);
	for (i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star(&xdrs,&sp);
	    for (j = 0; j < th.ndim; j++) {
		pos[j][i] = sp.pos[j];
		vel[j][i] = sp.vel[j];
		}
	    mass[i] = sp.mass;
	    eps[i] = sp.eps;
	    phi[i] = sp.phi;
	    metals[i] = sp.phi;
	    tform[i] = sp.tform;
	    }
	DBMkDir(dbfile,"star");
	DBSetDir(dbfile,"/star");
	DBPutPointmesh(dbfile,"star_pos",th.ndim,pos,th.nstar,DB_FLOAT,NULL);
	DBPutPointvar(dbfile,"star_vel","star_pos",th.ndim,vel,th.nstar,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"star_mass","star_pos",mass,th.nstar,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"star_eps","star_pos",eps,th.nstar,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"star_phi","star_pos",phi,th.nstar,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"star_metals","star_pos",metals,th.nstar,DB_FLOAT,NULL);
	DBPutPointvar1(dbfile,"star_tfrom","star_pos",tform,th.nstar,DB_FLOAT,NULL);
	DBPutDefvars(dbfile,"variables",nvar,star_names,types,star_defs,NULL);
	DBSetDir(dbfile,"/");
	}
    /*
    ** Clean up a bit
    */
    xdr_destroy(&xdrs);
    free(pos);
    free(vel);
    free(eps);
    free(phi);
    free(rho);
    free(temp);
    free(metals);
    free(tform);
    /*
    ** Now check if there are some additional arrays (use mass array for that)
    */
    if (denfile == 1) {
	file = fopen(denfilename,"r");
	assert(file != NULL);
	assert(fscanf(file,"%d",&i) == 1);
	assert(i == th.ntotal);
	if (th.ngas > 0) {
	    mass = realloc(mass,th.ngas*sizeof(float));
	    for (i = 0; i < th.ngas; i++) {
		assert(fscanf(file,"%f",&mass[i]) == 1);
		}
	    DBSetDir(dbfile,"/gas");
	    DBPutPointvar1(dbfile,"gas_den","gas_pos",mass,th.ngas,DB_FLOAT,NULL);
	    DBSetDir(dbfile,"/");
	    }
	if (th.ndark > 0) {
	    mass = realloc(mass,th.ndark*sizeof(float));
	    for (i = 0; i < th.ndark; i++) {
		assert(fscanf(file,"%f",&mass[i]) == 1);
		}
	    DBSetDir(dbfile,"/dark");
	    DBPutPointvar1(dbfile,"dark_den","dark_pos",mass,th.ndark,DB_FLOAT,NULL);
	    DBSetDir(dbfile,"/");
	    }
	if (th.nstar > 0) {
	    mass = realloc(mass,th.nstar*sizeof(float));
	    for (i = 0; i < th.nstar; i++) {
		assert(fscanf(file,"%f",&mass[i]) == 1);
		}
	    DBSetDir(dbfile,"/star");
	    DBPutPointvar1(dbfile,"star_den","star_pos",mass,th.nstar,DB_FLOAT,NULL);
	    DBSetDir(dbfile,"/");
	    }
	fclose(file);
	}
    if (psdfile == 1) {
	file = fopen(psdfilename,"r");
	assert(file != NULL);
	assert(fscanf(file,"%d",&i) == 1);
	assert(i == th.ntotal);
	if (th.ngas > 0) {
	    mass = realloc(mass,th.ngas*sizeof(float));
	    for (i = 0; i < th.ngas; i++) {
		assert(fscanf(file,"%f",&mass[i]) == 1);
		}
	    DBSetDir(dbfile,"/gas");
	    DBPutPointvar1(dbfile,"gas_psd","gas_pos",mass,th.ngas,DB_FLOAT,NULL);
	    DBSetDir(dbfile,"/");
	    }
	if (th.ndark > 0) {
	    mass = realloc(mass,th.ndark*sizeof(float));
	    for (i = 0; i < th.ndark; i++) {
		assert(fscanf(file,"%f",&mass[i]) == 1);
		}
	    DBSetDir(dbfile,"/dark");
	    DBPutPointvar1(dbfile,"dark_psd","dark_pos",mass,th.ndark,DB_FLOAT,NULL);
	    DBSetDir(dbfile,"/");
	    }
	if (th.nstar > 0) {
	    mass = realloc(mass,th.nstar*sizeof(float));
	    for (i = 0; i < th.nstar; i++) {
		assert(fscanf(file,"%f",&mass[i]) == 1);
		}
	    DBSetDir(dbfile,"/star");
	    DBPutPointvar1(dbfile,"star_psd","star_pos",mass,th.nstar,DB_FLOAT,NULL);
	    DBSetDir(dbfile,"/ ");
	    }
	fclose(file);
	}
    /*
    ** Finish up and write some output
    */
    DBClose(dbfile);
    free(mass);
    fprintf(stderr,"Time: %g Ntotal: %d Ngas: %d Ndark: %d Nstar: %d\n",
	    th.time,th.ntotal,th.ngas,th.ndark,th.nstar);
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program converts tipsy standard binary format to silo format\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Please specify the following parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-o <name>   : output file in silo format\n");
    fprintf(stderr,"-den <name> : density file in tipsy array format\n");
    fprintf(stderr,"-psd <name> : phase space density file in tipsy array format\n");
    fprintf(stderr,"< <name>    : input file in tipsy standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

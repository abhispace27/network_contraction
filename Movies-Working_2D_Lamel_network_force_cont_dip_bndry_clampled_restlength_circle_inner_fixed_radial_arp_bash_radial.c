// Conjugate gradient to relax an elastic network
// **********************************************
// this code applies a restlength change in bonds
// near the force dipole and then we
// let the energy minimizer do its magic.
// **********************************************

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "RAN.c"	//Random number generator
#include "nrutil_jen.h" //NR lib

//*************************************************
//**************** Preprocessor Stuff *************
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define sign(N,P) ((P) >= 0 ? (fabs(N)):(-fabs(N)))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define MOD(P,M) ((P%M)+M)%M
//#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);
#define M_PI 3.1415926535897932//the rest 3846264338327


//*************************************************
//************** Lattice parameters ***************
#define L 64		//L CAN ONLY BE EVEN THE SUBROUTINES WILL NOT WORK OTHERWISE
#define N (L*L)
#define mu 1.0	//This control the spring between n.n

//#define mu_c 1.0                                           // this is the compressive mu

//#define kappa1 1.0e-6	//This controls the ARP spring
//#define kappa1 0.0	//This controls the ARP spring
#define kappa1 100.0	//This controls the ARP spring
#define arp_lim M_PI/180.0  // this is the limit after which ar is activated

//#define kappa2 1.0e-6	//This controls the filament bending.
//#define kappa2 0.0	//This controls the filament bending.

#define THETA M_PI*(60.0/180.0)	//Define the Arp 2/3 angular bond interaction rest position for a perfect lattice.
//#define PBOND 0.64
//#define PBOND 0.61

#define PARP 1.0     // probablity of arp springs
#define PERT 0.01
#define PERT1 0


#define UNITC 3.464101615//3.0*sqrt(2.0)  // previous written by David. Mine: 2*sqrt(3)
//*************************************************
//************* Subroutine Parameters *************
#define NRANSI
#define TOL 1.0e-7	//This sets the tolerance for finding the minimum

//#define DSTEP 1.0e-6	//This sets the episilon in the derivative

#define ITMAXX 1000000	//Conj. Grad. Max Iterations 
#define ITMAX 1000000	//dbrent Max iterations

#define EPS 1.0e-10
#define ZEPS 1.0e-10
#define GOLD 1.618034	//Golden Ratio for mnbrak
#define CGOLD 0.3819660
#define GLIMIT 100.0
#define TINY 1.0e-20

#define EMPTY -(3*N)-1	//Percolate subroutine

// new vriables defined by abhinav feb 20
#define ext_flag 0 	// entensile flag: 1= extensile dipoles; 0= contractile dipoles 

//double Ff[num_node];   // stores the value of forces
// dist_val is the distance between nodes of a dipole
// xcomp and ycomp are the factors that determine dipole force components 
double dist_val, xdist_val, ydist_val, xcomp, ycomp;

// dynamic array initiation:
int* dip_node=NULL; 
double* Ff=NULL;
int* loc_center = NULL;
//int loc_center;

int* dip_node_unique = NULL;   // array of unique dipole nodes 
int num_dip_unique;

int num_center;  // number of isotropic contraction points
int num_dip, num_node;

//**************************************************
//**************************************************

long ranseed1, ranseed2;

int connect[N][6],occupation[6*N],arpoccupation[6*N],dip_list[2];

double strain_bonds[6*N];   // stores values of strains of all bonds

double x[N],y[N],pr[2*N],subpr[2*N];
int filenum, outfilenum, outfilenum1, outfilenum2, tstep;

double brakbx=2.0*PERT,brakax=0.0;      // do not change unless really needed
//double brakbx=100.0,brakax=0.0;    // edit by abhinav; trying different values

int func_flag = 0;


int en_flag = 0;  // used to store energy values from conjugate method so that i can generate animation
int outfilenum3, outfilenum4; 	// used to store energy values from conjugate method so that i can generate animation
int outfilenum5;  // used to number 
int outfilenum_strain = 0;
int movie_flag = 1;  // to create movie of all conjugate gradient movements

int coord_num[N];  // cooordination number of the nodes

double ES[6*N];    // global so that we can write to file at the end

// making energies global so we can write them to file from wherever in the code.
double Eelastic,Ebend,Filbend,ETOT, Eforce;
double Etension, Etension_dip, Ecompression, Ebuckle;

int outfilenum_cg_en = 0;
int outfilenum_cg = 0;    // for file name containing location of nodes per cg step
int outfilenum_force = 0;

int times;   // making it global because it will go in filenames
double dip_move;

//double rlen_max = 0.8;
double rlen_max = 0.9;
double rlen;
double delta_bond_len = 0;

double fxi[2*N];   // this array stores forces; created because in this array boundary forces are not zero.
double rlen_arr[6*N];   // stores rest lengths of all bonds
int rlen_outfilenum = 0;   // to update rlen output file name
int en_node_outfilenum = 0;  // for energy per node filename

int outer_node[N];  // stores 0 or 1; 0: for nodes inside the circular boundary; 1: for the nodes on the boundary outside it

// the following are variables that change the location of the dipoles without changing the network
int srand_flag = 0;
// int srand_arg = 114;
int srand_arg;

// following are variables that change the outer network while the inner network is kept same
// 111,111   111,112     667720
//int diff_network1 = 111111;   // using the same variables as for changing dipole locations because the folders are already created

//int diff_network2 = 111111;


//****************************************************************
//******************** FOR BASH SCRIPT ***************************
//double mu_c_bash, PBOND_bash, kappa2_bash;
double mu_c, PBOND, kappa2;
int num_center_bash;
char arg_txt[255];
int diff_network1, diff_network2;

//****************************************************************
// test for writing the forces to file
double FILAx[N*3],FILAy[N*3];


//****************************************************************
//********************OUTPUT FILES BELOW**************************
//****************************************************************
//
// Main idea:
// 1. make_output_filename() builds the full filename.
// 2. open_output_file() safely opens the file.
// 3. Each output_* function only writes its own data.
//
//****************************************************************

#define OUT_ROOT "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial"
#define FILENAME_SIZE 1000

// Example:
// 1.0e-6  -> e-6
// 2.0e-6  -> 2e-6
// 5.0e-7  -> 5e-7
// 0       -> 0
//
void format_kappa2_folder(double val, char *buf, size_t n)
{
    if (val == 0.0) {
        snprintf(buf, n, "0");
    }
    else if (val == 1.0e-7) {
        snprintf(buf, n, "e-7");
    }
    else if (val == 1.0e-6) {
        snprintf(buf, n, "e-6");
    }
    else if (val == 1.0e-5) {
        snprintf(buf, n, "e-5");
    }
    else if (val == 1.0e-4) {
        snprintf(buf, n, "e-4");
    }
    else if (val == 1.0e-3) {
        snprintf(buf, n, "e-3");
	}
    else if (val == 1.0e-2) {
        snprintf(buf, n, "e-2");
    }
    else {
        snprintf(buf, n, "%.0e", val);
    }
}

// Adds L prefix used by some old filenames.
void format_L_prefix(char *buf, size_t n)
{
    if (L == 96) {
        snprintf(buf, n, "96_");
    }
    else if (L == 128) {
        snprintf(buf, n, "128_");
    }
    else {
        snprintf(buf, n, "");
    }
}


// Build the folder before the actual filename.
//
// For PBOND == 1:
// OUT_ROOT/subdir/
//
// For PBOND < 1:
// OUT_ROOT/kappa2_x/ranseed1,ranseed2/srand_arg/diff_network1,diff_network2/subdir/

void make_output_dir(char *dir, size_t n, const char *subdir)
{
    char kappa_txt[64];

    if (PBOND == 1.0) {
        snprintf(dir, n, "%s/%s", OUT_ROOT, subdir);
    }
    else {
        format_kappa2_folder(kappa2, kappa_txt, sizeof(kappa_txt));

        snprintf(dir, n,
                 "%s/kappa2_%s/%ld,%ld/%d/%d,%d/%s",
                 OUT_ROOT,
                 kappa_txt,
                 ranseed1, ranseed2,
                 srand_arg,
                 diff_network1, diff_network2,
                 subdir);
    }
}


// Builds the common parameter part of the filename.
void make_common_suffix(char *suffix, size_t n, int counter, int use_counter)
{
    char L_txt[32];

    format_L_prefix(L_txt, sizeof(L_txt));

    if (use_counter) {
        snprintf(suffix, n,
                 "%s%d_%d_%0.2f_%0.2e_%0.2e__%0.4f_%0.4f_%0.4f_%d_%d",
                 L_txt,
                 num_center,
                 num_dip,
                 PBOND,
                 TOL,
                 kappa2,
                 rlen_max,
                 mu,
                 mu_c,
                 counter,
                 times);
    }
    else {
        snprintf(suffix, n,
                 "%s%d_%d_%0.2f_%0.2e_%0.2e__%0.4f_%0.4f_%0.4f_%d",
                 L_txt,
                 num_center,
                 num_dip,
                 PBOND,
                 TOL,
                 kappa2,
                 rlen_max,
                 mu,
                 mu_c,
                 times);
    }
}


// General filename builder.
void make_output_filename(
    char *filename,
    size_t n,
    const char *subdir,
    const char *prefix,
    int counter,
    int use_counter,
    const char *extension
)
{
    char dir[FILENAME_SIZE];
    char suffix[FILENAME_SIZE];

    make_output_dir(dir, sizeof(dir), subdir);
    make_common_suffix(suffix, sizeof(suffix), counter, use_counter);

    snprintf(filename, n, "%s/%s%s%s", dir, prefix, suffix, extension);
}


// Special filename builder for node energy.
void make_energy_node_filename(char *filename, size_t n, int counter)
{
    char dir[FILENAME_SIZE];
    char suffix[FILENAME_SIZE];
    char L_txt[32];

    format_L_prefix(L_txt, sizeof(L_txt));

    if (PBOND == 1.0) {
        snprintf(dir, sizeof(dir), "%s/energy/node", OUT_ROOT);

        snprintf(suffix, sizeof(suffix),
                 "%ssrand_%d_%d_%d_%0.2f_%0.2e_%0.2e_%0.4f_%0.4f_%0.4f_%d_%d",
                 L_txt,
                 srand_arg,
                 num_center,
                 num_dip,
                 PBOND,
                 TOL,
                 kappa2,
                 rlen_max,
                 mu,
                 mu_c,
                 counter,
                 times);

        snprintf(filename, n, "%s/Lattice_node_%s_force.txt", dir, suffix);
    }
    else {
        make_output_filename(filename, n,
                             "energy/node",
                             "Lattice_node_",
                             counter,
                             1,
                             "_force.txt");
    }
}


// Safe fopen check.
FILE *open_output_file(const char *filename, const char *mode)
{
    FILE *fout = fopen(filename, mode);

    if (fout == NULL) {
        fprintf(stderr, "\nERROR: could not open file:\n%s\n", filename);
        exit(EXIT_FAILURE);
    }

    return fout;
}


//****************************************************************
// Positions + occupation + coordination number
//****************************************************************

void output_sub(double pos[])
{
    int i, j;
    char filename[FILENAME_SIZE];

    make_output_filename(filename, sizeof(filename),
                         "txt/strain",
                         "Lattice_",
                         outfilenum,
                         1,
                         "_force.txt");

    outfilenum++;

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < N; i++) {
        fprintf(fout, "%d %20.15f %20.15f ", i, pos[i], pos[i + N]);

        for (j = 0; j < 6; j++) {
            fprintf(fout, " %d", occupation[i * 6 + j]);
        }

        fprintf(fout, " %d\n", coord_num[i]);
    }

    fclose(fout);

    printf("\nFilename is: %s\n", filename);
}


//****************************************************************
// Forces on nodes
//****************************************************************

void force_output(double force[])
{
    int i;
    char filename[FILENAME_SIZE];

    make_output_filename(filename, sizeof(filename),
                         "txt/force",
                         "Force_",
                         outfilenum_force,
                         1,
                         "_force.txt");

    outfilenum_force++;

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < N; i++) {
        fprintf(fout, "%15.5e \t %15.5e\n", force[i], force[i + N]);
    }

    fclose(fout);
}


//****************************************************************
// Bending forces
//****************************************************************

void bend_force_output()
{
    int i, j;
    char filename[FILENAME_SIZE];

    make_output_filename(filename, sizeof(filename),
                         "txt/area",
                         "bend_force_",
                         0,
                         0,
                         "_force.txt");

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < N; i++) {
        fprintf(fout, "%12d \t", i);

        for (j = 0; j < 3; j++) {
            fprintf(fout, "%15.5e \t %15.5e \t",
                    FILAx[3 * i + j],
                    FILAy[3 * i + j]);
        }

        fprintf(fout, "\n");
    }

    fclose(fout);
}


//****************************************************************
// Rest lengths
//****************************************************************

void rlen_output()
{
    int i, j;
    char filename[FILENAME_SIZE];

    make_output_filename(filename, sizeof(filename),
                         "txt/rlen",
                         "rlen_",
                         rlen_outfilenum,
                         1,
                         ".txt");

    rlen_outfilenum++;

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < N; i++) {
        for (j = 0; j < 6; j++) {
            fprintf(fout, "%8.4f ", rlen_arr[i * 6 + j]);
        }
        fprintf(fout, "\n");
    }

    fclose(fout);
}


//****************************************************************
// Generic node-list writer
// Used by dipole, unique dipole, boundary, and inner node outputs
//****************************************************************

void output_node_list(
    const char *subdir,
    const char *prefix,
    const int nodes[],
    int count
)
{
    int i;
    char filename[FILENAME_SIZE];

    make_output_filename(filename, sizeof(filename),
                         subdir,
                         prefix,
                         0,
                         0,
                         "_force.txt");

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < count; i++) {
        fprintf(fout, " %d %d \n", i, nodes[i]);
    }

    fclose(fout);
}


//****************************************************************
//****************************************************************

void output_dip_nodes()
{
    output_node_list("txt/dip_nodes",
                     "Dipole_nodes_",
                     dip_node,
                     num_node);
}


//****************************************************************
//****************************************************************

void output_unique_dip_nodes()
{
    output_node_list("txt/dip_nodes",
                     "Unique_dipole_nodes_",
                     dip_node_unique,
                     num_dip_unique);
}


//****************************************************************
//****************************************************************

void output_boundary_nodes(int boundary_nodes[], int temp_count)
{
    output_node_list("txt/area",
                     "Boundary_nodes_",
                     boundary_nodes,
                     temp_count);
}


//****************************************************************
//****************************************************************

void output_inner_nodes(int inner_nodes[], int inner_count)
{
    output_node_list("txt/area",
                     "Inner_nodes_",
                     inner_nodes,
                     inner_count);
}


//****************************************************************
// Outer nodes
//****************************************************************

void output_outer_nodes()
{
    int i;
    char filename[FILENAME_SIZE];

    make_output_filename(filename, sizeof(filename),
                         "txt/area",
                         "Outer_nodes_",
                         0,
                         0,
                         "_force.txt");

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < N; i++) {
        if (outer_node[i] == 1) {
            fprintf(fout, " %6d \n", i);
        }
    }

    fclose(fout);
}


//****************************************************************
// Connectivity output
//****************************************************************

void output_connect(double pos[])
{
    int i, j;
    int dummy;
    char filename[FILENAME_SIZE];

    make_output_filename(filename, sizeof(filename),
                         "txt/strain",
                         "Lattice_connect_",
                         outfilenum5,
                         1,
                         "_force.txt");

    outfilenum5++;

    printf("\nfilename: %s\n", filename);

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < N; i++) {
        fprintf(fout, "%5d %014.10f %014.10f ", i, pos[i], pos[i + N]);

        for (j = 0; j < 6; j++) {
            if (occupation[connect[i][j] * 6 + (j + 3) % 6] == 1) {
                fprintf(fout, " %5d", connect[i][j]);
            }
            else {
                if (L == 64) {
                    dummy = 9999;
                }
                else {
                    dummy = -1;
                }

                fprintf(fout, " %5d", dummy);
            }
        }

        fprintf(fout, " %d\n", coord_num[i]);
    }

    fclose(fout);
}


//****************************************************************
// Strain and stretching energy per bond
//****************************************************************

void output_strain_sub()
{
    int i, j;
    char filename[FILENAME_SIZE];

    make_output_filename(filename, sizeof(filename),
                         "txt/strain",
                         "Strain_Lattice_",
                         outfilenum_strain,
                         1,
                         "_force.txt");

    outfilenum_strain++;

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < N; i++) {
        fprintf(fout, "%d ", i);

        for (j = 0; j < 6; j++) {
            fprintf(fout, " %15.10f", strain_bonds[i * 6 + j]);
        }

        for (j = 0; j < 6; j++) {
            fprintf(fout, " %15.10e", ES[i * 6 + j]);
        }

        fprintf(fout, "\n");
    }

    fclose(fout);
}


//****************************************************************
// Total system energy
//****************************************************************

void output_energy_sub(
    double en_tot,
    double en_elastic,
    double en_comp,
    double en_bend,
    double en_arp
)
{
    char filename[FILENAME_SIZE];
    const char *mode;

    make_output_filename(filename, sizeof(filename),
                         "energy/strain",
                         "Lattice_",
                         0,
                         0,
                         "_force.txt");

    if (outfilenum1 == 0) {
        mode = "w";
    }
    else {
        mode = "a";
    }

    outfilenum1++;

    FILE *fout = open_output_file(filename, mode);

    fprintf(fout, "%15.5e %15.5e %15.5e %15.5e %15.5e\n",
            en_tot,
            en_elastic,
            en_comp,
            en_bend,
            en_arp);

    fclose(fout);
}


//****************************************************************
// Energy at each node
//****************************************************************

void output_energy_node_sub(double ES[], double FILB[])
{
    int i, j;
    char filename[FILENAME_SIZE];

    make_energy_node_filename(filename, sizeof(filename), en_node_outfilenum);

    en_node_outfilenum++;

    FILE *fout = open_output_file(filename, "w");

    for (i = 0; i < N; i++) {
        fprintf(fout, "%10d ", i);

        // stretching energy
        for (j = 0; j < 6; j++) {
            fprintf(fout, "%15.5e ", ES[6 * i + j]);
        }

        // bending energy
        for (j = 0; j < 3; j++) {
            fprintf(fout, "%15.5e ", FILB[3 * i + j]);
        }

        fprintf(fout, "\n");
    }

    fclose(fout);
}


//****************************************************************
// Creates output filename for storing bend/angle data
//
// Important:
// This returns malloc'ed memory.
//****************************************************************

char *output_bend_name()
{
    char *filename = malloc(sizeof(char) * FILENAME_SIZE);

    if (filename == NULL) {
        fprintf(stderr, "ERROR: malloc failed in output_bend_name()\n");
        exit(EXIT_FAILURE);
    }

    make_output_filename(filename, FILENAME_SIZE,
                         "txt/bend/strain",
                         "Lattice_bend_",
                         outfilenum2,
                         1,
                         "_force.txt");

    outfilenum2++;

    return filename;
}

//****************************************************************
//********************OUTPUT FILES ABOVE**************************

//****************************************************************
//****************************************************************

void eraserhead(){
	int i,j;

	for(i=0;i<N;i++){
		for(j=0;j<6;j++){
			occupation[i*6+j]=0;
		}
	}
	
}

double dangle(){

int i,j,k,p,temp1,hold1,hold2,tick;
double count;


i=j=p=temp1=hold1=hold2=0;
tick=1;

	while(tick>0){
//		dlattice(pr);		
		tick=0;
		for(i=0;i<N;i++){
			for(j=0;j<6;j++){
				if(occupation[i*6+j]==1 && temp1<2){
					temp1++;					
					hold1=i*6+j;
					hold2=connect[i][j]*6+(j+3)%6;						
										
				}
			}
			if(temp1==1){
				occupation[hold1]=occupation[hold2]=0;
				
				temp1=0;
				tick++;
//				printf("%d\n",tick);	
			}	
			if(temp1>1){
				temp1=0;
							
			}			
		}
	}

for(i=0;i<N;i++){
		for(j=0;j<6;j++){
			count=count+occupation[i*6+j];
		}
	}

tick=0;
k=0;
	for(i=0;i<N;i++){
		k++;
		for(j=0;j<6;j++){
			if(i<N-L && j>3 && j<5){
				if(occupation[i*6+j]==0){
					if( occupation[i*6+(j+1)%6]==0){
						tick++;
						//printf("%d %d\n",k/L,tick);
					}
				}
				
				if(k%L==0){
					if(tick/L==1){
//						printf("ticky-tick\n");	
//						printf("%d %d\n",k/L,tick);
//						dlattice(pr);
//						eraserhead();						
						tick=0;
						
					}
					tick=0;
				}
			}
		}	
	}
	return count/(6.0*N-4.0*L);
}


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&& CONJUGATE GRADIENT &&&&&&&&&&&&&&&&&&&&&&&&&

// This subroutine attempts to find the global minimum of a function
// "func()" that is specified by the user. To find the global minimum
// it assumes that it exists and can be found by minimizing func along 
// each conjugate direction in coordinate space of size 2*N. It returns 
// the final configuration of all the coordinates once the minimum is 
// within some tolerance TOL. 

void frprmn( double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), void (*dfunc)(double [], double []) )
{	// p[] - is the lattice configuration which is 2*number of unconstrained degrees of freedom
	// q[] is the same configuration with the boundaries distorted
	// c is the connection array(matrix)
	// o is the bond occupation array(matrix)
	// n is the size of the lattice N - number of fixed sites 
	// ftol the tolernce for the final minimum when found
	// *iter is the number of iterations that the algorithm need to find the minimum (good for bench marking)
	// *fret is the returned value of the minimum found by dbrent after the mine minimazation
	 

	void dlinmin( double p[], double xi[], int n, double *fret, double (*func)(double [] ), void (*dfunc)(double [], double []) );

	int j,its;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;
//	printf("Hello, you made it in! Conj Grad!\n");
//	printf("\n");
	g=dvector(1,n);
	h=dvector(1,n);
	xi=dvector(1,n);
	fp=(*func)(p);
//	printf("IN THE BEGINING func = %f\n",fp);
	(*dfunc)(p,xi);

// for finding minima and maxima of force values to study extremely large values in the enrgies
	int ii;
	double min, max;  

	for (j=1;j<=n;j++) { 
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
//		printf(" step j: %d xi: %f g: %f h: %f p[j]: %f \n", j, xi[j], g[j], h[j], p[j]);
	}
	for (its=1;its<=ITMAXX;its++) {
		*iter=its;

//		printf("\n BEGIN ----- before dlinmin fp: %20.15f fret: %20.15f \n", fp, *fret);
		dlinmin(p,xi,n,fret,func,dfunc); //this use derivative information to find the min
//		printf("\n 2.0*fabs(*fret-fp): %20.15f \n", 2.0*fabs(*fret-fp));
//		printf("\n ftol*(fabs(*fret)+fabs(fp)+EPS): %20.15f \n", ftol*(fabs(*fret)+fabs(fp)+EPS));
//		printf("\n END   ----- after dlinmin fp: %20.15f fret: %20.15f \n", fp, *fret);



		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) { //check if we got it!
//			for (j=1;j<n; j++)
//			{
//				printf("\n xi[%d]: %0.10f ", j, xi[j]);				/* code */
//			}

		min=max=xi[0];
	    for(ii=1; ii<n; ii++)
	    {
    	    if(min>xi[ii])
				min=xi[ii];   
		   	if(max<xi[ii])
		    	max=xi[ii];       
	    }

    	printf("minimum of array is : %e", min);
        printf("\nmaximum of array is : %e", max);
		
		printf("\n n: %d \n", n);

			printf("\n 1. Found the minima:  %f after:  %d steps \n", fp, its);

			force_output(fxi);

			free_dvector(xi,1,n);
			free_dvector(h,1,n);
			free_dvector(g,1,n);


//			FREEALL
			return;
		}
//if not then set the new value of fp to that of the last value of what dlinmin spits out
		fp=*fret;
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {// build a new conjugate direction to test // changed from j=1 j<=n
			gg += g[j]*g[j];
			//dgg +=xi[j]*xi[j];	
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (dgg == 0.0) { //changed from gg to dgg
			free_dvector(xi,1,n);
			free_dvector(h,1,n);
			free_dvector(g,1,n);

			printf("\n 2. Found the minima:  %f after:  %d steps \n", fp, its);

//			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
//			printf(" step j: %d xi: %f g: %f h: %f gam: %f \n", j, xi[j], g[j], h[j], gam);		
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
//			printf(" step j: %d xi: %f g: %f h: %f gam: %f , p[j]: %f \n", j, xi[j], g[j], h[j], gam, p[j]);		
		}
//		printf("One more time!\n");
	
	}
	printf("\n in frprmn line 1797. not a good place to be \n");

	nrerror("Too many iterations in frprmn");
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&& LINE MINIMAZATION &&&&&&&&&&&&&&&&&&&&

//This sub routine is used by conjugate gradient to find the minimum along one direction
//in the parameter space. Which in this case is 2*N dimensions.


int ncom;
double *pcom,*xicom;

double (*nrfunc)( double []);
void (*nrdfun)( double [], double []);

void dlinmin( double p[], double xi[], int n, double *fret, double (*func)(double []), void (*dfunc)(double [], double []) )

//void dlinmin( double p[], double xi[], int n, double *fret, double (*func)(double []), void (*dfunc)(double [], double []) )
{
	double dbrent( double ax, double bx, double cx, double(*f)(double ), double (*df)(double), double tol, double *xmin );

	double f1dim(double x);
	
	double df1dim(double x);

	void mnbrak( double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double) );
	
	int j;
	double cx,bx,ax;
	double xmin,fx,fb,fa;

//	printf("We are in dlinmin\n");

	ncom=n;
	pcom=dvector(1,n);
	xicom=dvector(1,n);
	nrfunc=func;
	nrdfun=dfunc;
	for(j=1;j<=n;j++)//changed from j=1 and j<=n
		{
			pcom[j]=p[j-1]; // j-1 used to be j, need for function definition. NR requires everything start at j=1
			xicom[j]=xi[j];
//			printf(" step j: %d pcom: %f xicom: %f \n", j, pcom[j], xicom[j]);
		}
	ax=brakax;			//Set bracket size to previous max -- this is the guess for the minima bracketing
	bx=brakbx;	
	cx=0.0;
	
//	printf(" before mnbrak ax: %f bx: %f cx: %f \n", ax, bx, cx);

	mnbrak(&ax,&bx,&cx,&fa,&fx,&fb,f1dim);

//	printf(" after mnbrak ax: %f bx: %f cx: %f \n", ax, bx, cx);
	
//	printf("\n before dbrent fret: %f \n", *fret);

	*fret=dbrent(ax,bx,cx,f1dim,df1dim,TOL,&xmin);

// 	printf("\n after dbrent fret: %f \n", *fret);


// 	printf(" ax: %f bx: %f cx: %f \n", ax, bx, cx);

 	brakbx=bx;		//Store previous bracket max
 	brakax=ax;		

 	for(j=1;j<=n;j++)	//Construct the dvector to return //changed from j=1 j<=n
 	{
 		xi[j] *= xmin;
 		p[j-1] += xi[j]; //need for function definition. NR requires everything start at j=1, but func starts at j=0.

 //		printf(" step j: %d xi: %f p: %f \n", j, xi[j], p[j]);
 	}
	
 	free_dvector(xicom,1,n);
 	free_dvector(pcom,1,n);
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&& BRACKETING SUBROUTINE &&&&&&&&&&&&&&&&&&&&&&&

// This subroutine loactes the minimum within a bracket of values of
// the function func() and the absissa of the line dlinmin is attempting
// to minimize along. It returns 3 valuse of func() and three points ax,bx,cx
// where the minimum is approximatley loacted using a parabolic fit. 

void mnbrak( double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double) ) 
{
	double ulim,u,r,q,fu,dum;
	double t,tt,tiny;
	*fa=(*func)(*ax); //this is actually calling f1dim which does take a single number
	*fb=(*func)(*bx);
//	printf("\n fa: %f fb: %f fc: %f ax: %f bx: %f cx: %f dum: %f \n", *fa, *fb, *fc, *ax, *bx, *cx, dum);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
//	printf("\n fa: %f fb: %f fc: %f ax: %f bx: %f cx: %f dum: %f \n", *fa, *fb, *fc, *ax, *bx, *cx, dum);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		tt=fabs(q-r);
		tiny=TINY;		
		t=fmax(tt,tiny);	
	
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*sign(t,q-r));

		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
//		printf("\n fa: %f fb: %f fc: %f ax: %f bx: %f cx: %f dum: %f \n", *fa, *fb, *fc, *ax, *bx, *cx, dum);
	}
}
	


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&& 1D MINIMIZATION SUBROUTINE &&&&&&&&&&&&&&&&&&

// dbrent uses the information that is collected by mbrack and finds the
// minimum of func along the line within the bracket by evaluating func and dfunc
// by calling f1dim and df1dim. 

double dbrent( double ax, double bx, double cx, double (*f)(double),double (*df)(double), double tol, double *xmin )
{
	int iter,ok1,ok2;
	double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);//this calls f1dim
	dw=dv=dx=(*df)(x);//this calls df1dim
//	printf("Now we made it to dbrent func = %f\n",fw);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx) d1=(w-x)*dx/(dx-dw);
			if (dv != dx) d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=sign(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu=(*f)(u);
		} else {
			u=x+sign(tol1,d);
//printf("inside dbrent calling f1dim again\n");
			fu=(*f)(u);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
//printf("inside dbrent calling df1dim again\n");
		du=(*df)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
	nrerror("Too many iterations in routine dbrent");
	return 0.0;
}


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&& 1D DERIVATIVE &&&&&&&&&&&&&&&&&&&&&&&

// This subroutine which is called by dbrent finds the derivative
// using dfunc of func when it is being minimized along one direction

double df1dim(double x)
{
	int j;
	double df1=0.0;
//	double xt[ncom],*df;
	double *xt,*df;
//printf("%d\n",ncom);
	xt=dvector(0,ncom);
	df=dvector(1,ncom);
	
	for(j=1;j<=ncom;j++){
		 xt[j-1]=pcom[j]+x*xicom[j]; // j-1 needed for function definition. NR requires everything start at j=1
	}

	(*nrdfun)(xt,df);
	
	for(j=1;j<=ncom;j++) df1 +=df[j]*xicom[j];

	free_dvector(df,1,ncom);
	free_dvector(xt,0,ncom);
	return df1;
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&& 1D FUNC EVALUATION &&&&&&&&&&&&&&&&&&

// This routine evaluates func along the direction that dbrent
// is trying to locate the minimum along. It's value is used to
// bracket the minimum. 


double f1dim(double x) // this routine talks to mnbrak and dbrent. It evalutates func along one direction, instead of calling func itself
{
	int j;
	double f;
//	double xt[ncom];
	double *xt;
	xt=dvector(0,ncom);			
	for (j=1;j<=ncom;j++) 
	{	//take the initial point p[] and move by an amount x in the direction xicom[]
		xt[j-1]=pcom[j]+x*xicom[j];  // j-1 needed for function definition. NR requires everything start at j=1
//		printf("pcom=%f, xicom=%f\n",pcom[j],xicom[j]);
	}
	f=(*nrfunc)(xt);	//this calls func.
//printf("Were in f1dim and func = %f\n",f);
	free_dvector(xt,0,ncom);
	return f;
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&& FUNC &&&&&&&&&&&&&&&&&&&&&&&&&&&&

//This subroutine calculates the energy of the lattice given is 
//configuration which is stored in x[N] and y[N]. Returns a single
//value.


double func(double pit[])//,
{
	int k,i,j,m,bone,btwo,bthree,arp;
//	double Eelastic,Ebend,Filbend,ETOT, Eforce;
//	double Etension, Etension_dip, Ecompression, Ebuckle;
	int dip_true, mu_flag;    // flag to check if we are at a dipole node to discriminate between buckling and compression
	double deltaX,deltaY;
	double EB[6*N];
//	double ES[6*N],FILB[3*N];   // commenting this to make ES a global variable to write to file
	double FILB[3*N];
	double newx[N],newy[N];
	double delxij,delxik,delyij,delyik,delxjk,delyjk;
	double THETAijk,THETAjik,THETAjki,DOTijk,DOTjik,DOTikj,CROSSjik;
	double Rij,Rik,Rjk;

	double THETAjik_half, sin_THETAjik_half;

	double Ecomp;

	int d_i;

	int jj;

	int i_prev, i_new;  // used to write new line in bend file

	int dip_flag = 0;
	char* filename;
	FILE * fout;
	FILE * fout_en;
	int row;
	
	
	// Eelastic  is the energy stored in the lattice due to two-body interactions.
	// Ebend is the energy stored in the lattice due to the 3-body interactions.
	// deltaRx and deltaRy is the change in the bond lengths in the x and y directions respectively
	// F[3*N] is the energy of each of the bonds 3*N per site
	//delxij,delxik,delyij and delyik are the distances x and y between sites i,j and k. (3-body)
	//thetanew[] is the angle between neighboring bonds
	//BE[15*N] is the number 3-body interactions that each site [i] has 


/*	   [k]   *
	     \   |	       
	      *  *
	       \ |		
     [j] __ *__ [i] __  __ *
		 | \
		 *  *
		 |   \
		 *    *   */

	//Rij, Rik, Rjk are the bond lengths of the ij, ik and jk bonds respectively 

	
	Eelastic=0.0;
	Ebend=0.0;
	ETOT=0.0;
	Filbend=0.0;

	Etension=0.0;
	Etension_dip=0.0;
	Ecompression=0.0;
	Ebuckle=0.0;
	Ecomp=0.0;

// energy due to the force dipoles at nodes
	Eforce = 0.0;

//	j=0;

// condition below to keep ALL boundaries fixed.
	for(i=0;i<N;i++){
		
		dip_true = 0;    // flag for checking if we are at a dipole node

		if (outer_node[i] == 1){      // outer boundary
			newx[i]=pr[i];
			newy[i]=pr[i+N];
//			j++;
		}
		else{
			newx[i]=pit[i];
			newy[i]=pit[i+N];
//			j++;
		}

	}

	for(i=0;i<N;i++){
		for(j=0;j<6;j++){
		
			EB[i*6+j]=ES[i*6+j]=0.0;
			if(j/3==0){
				FILB[3*i+j]=0.0;
			}
		}
	}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!! BOND SPRING !!!!!!!!!!!!!!!!!!!!!!!!!!!!
if( mu > 0.0){
//printf("Not Supposed to be in here\n");
Eelastic=0.0;
Etension=0.0;
Etension_dip=0.0;
Ecompression=0.0;
Ebuckle=0.0;
	for ( i=0; i<N; i++ ) { 	
		for( j=0; j<6; j++ ) {
//			printf("\n In Eelastic loop i: %d j: %d \n", i, j);
			deltaX=deltaY=0.000000;
			Rij=0.000000;
			ES[6*i+j]=0.0000000;

			dip_true = 0;
			if(occupation[i*6+j]==1){

				deltaX=	newx[i]-(newx[ connect[i][j] ]);
//				deltaX=deltaX-( rint(deltaX/L)*L );
		
				deltaY=	newy[i]-(newy[ connect[i][j] ]);
//				deltaY=deltaY-( rint(deltaY/(L*sin(THETA)) )*L*sin(THETA) );
				Rij=sqrt( pow(deltaX,2)+pow(deltaY,2) );	

				if (func_flag == 1) strain_bonds[6*i+j] = Rij;   // storing strain to output

//				ES[6*i+j]=mu/4*pow( (Rij-1.0) ,2);

//				delta_bond_len = Rij-1.0;
				delta_bond_len = Rij-rlen_arr[i*6+j];

				ES[6*i+j]=mu/4*pow( delta_bond_len ,2);

			}

			if (delta_bond_len >= 0) Eelastic=Eelastic+ES[6*i+j];
			else Ecomp=Ecomp+ES[6*i+j];
			//Eelastic=Eelastic+ES[6*i+j];
				
		}
		
	}
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!! ARP ANGULAR SPRING !!!!!!!!!!!!!!!!!!!!!!

delxij=delyij=delxik=delyik=0.0;
if(kappa1 > 0.0){
	for( i=0;i<N;i++ ){
		for( j=0;j<6;j++ ){
						
			bone=occupation[i*6+j]; 
			btwo=occupation[i*6+(j+1)%6];  
			arp=arpoccupation[i*6+j]==1;			
//			if( bone==1 && btwo==1 && arp==1 ){// && j>0 && j<3 ){ //ANISO-exclude all spings except the one between 1 and 2
			if( bone==1 && btwo==1 ){ // removed arp check from above because arp is everywhere
//			
//			&&&&&&&&&&&&&&&& 1st Bond &&&&&&&&&&&&&&&&&				  	

				delxij = newx[i]-(newx[ connect[i][j] ]);	
//				delxij = delxij-( rint(delxij/L)*L );

				delyij = newy[i]-(newy[ connect[i][j] ]);
//				delyij = delyij-( rint(delyij/(L*sin(THETA)) )*L*sin(THETA) );

				Rij = sqrt( pow(delxij,2)+pow(delyij,2) );
//				printf("Rij=%f for sites %d and %d\n",Rij,i,connect[i][j]);
	
//			&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		
//			&&&&&&&&&&&&&&&& 2nd Bond &&&&&&&&&&&&&&&&&

				delxik = newx[i]-( newx[ connect[i][(j+1)%6] ]);
//				delxik = delxik-( rint(delxik/L)*L );

				delyik = newy[i]-( newy[ connect[i][(j+1)%6] ] );
//				delyik = delyik-( rint(delyik/(L*sin(THETA)) )*L*sin(THETA) );

				Rik = sqrt( pow(delxik,2)+pow(delyik,2) );
//				printf("Rik=%f for sites %d and %d\n",Rik,i,connect[i][(j+1)%6]);



//			&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	
//			&&&&&&&&&&&&&&&&& Compute the angle Between jik &&&&&&&&&&&&&&&			

//			printf("for sites %d,%d and sites %d,%d\n",i,connect[i][j],i,connect[i][(j+1)%6]);
/*a*/		
//DOT PROD
			THETAjik=acos((delxij*delxik+delyij*delyik)/(Rij*Rik));// Dot Product

			//EB[6*i+j]=kappa1/2*pow( (THETAjik-THETA),2 );
			//Ebend=Ebend+EB[6*i+j];

			if (THETAjik < arp_lim){    // checking if two bonds are about to cross beyond 1 degree separation

	//			printf("THETAjik=%f\n",THETAjik*180/M_PI);
	//			printf("\n");
	//			Ebend=Ebend+kappa1/2*pow( (THETAijk-THETA),2 );
				EB[6*i+j]=kappa1/2*pow( (THETAjik-arp_lim),2 );
				Ebend=Ebend+EB[6*i+j];

			}
			else {
				EB[6*i+j]=0;
				Ebend=Ebend+EB[6*i+j];
			}

				
			}

		}	
	}	
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!! FILAMENT ANGULAR SPRING !!!!!!!!!!!!!!!!!!!!!!

// opening bend file
if(func_flag == 1){
	printf("\n we here\n");
	filename = output_bend_name();
	printf("%s\n",filename);
	fout=fopen( filename,"w");
//	fd=fopen("Lattice.ps","w");	//open the file or create one if it doesn't exist, "w" is the 'write to' flag
}

if(kappa2 > 0.0){
Filbend=0.0;
	for( i=0;i<N;i++ ){
		for( j=0;j<3;j++ ){

			if (i == 0){
				i_prev = 0;             // used to write new line in bend file
			}

			CROSSjik=0.00000;
			THETAjik=0.00000;
			delxij=delyij=delxik=delyik=0.00000;
			Rij=Rik=0.000000;

			bone=occupation[i*6+j]; 
			btwo=occupation[i*6+(j+3)%6];  
//			if(  (i<L && connect[i][j]>=N-L)  || (i>=N-L &&  connect[i][(j+3)%6] < L) ){
//				FILB[3*i+j]=0.0;
//			}else 
			if( (bone==1 && btwo==1) ){  

//				&&&&&&&&&&&&&&&& 1st Bond &&&&&&&&&&&&&&&&&				  	

				delxij = newx[i]-(newx[ connect[i][j] ]);	
//				delxij = delxij-( rint(delxij/L)*L );

				delyij = newy[i]-(newy[ connect[i][j] ]);
//				delyij = delyij-( rint(delyij/(L*sin(THETA)) )*L*sin(THETA) );

				Rij = sqrt( pow(delxij,2)+pow(delyij,2) );

//				printf("Rij=%f for sites %d and %d\n",Rij,i,connect[i][j]);
//				printf("Bond=%d\n",j);
//				&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		
//				&&&&&&&&&&&&&&&& 2nd Bond &&&&&&&&&&&&&&&&&

				delxik = newx[i]-( newx[ connect[i][(j+3)%6] ] );
//				delxik = delxik-( rint(delxik/L)*L );

				delyik = newy[i]-( newy[ connect[i][(j+3)%6] ] );
//				delyik = delyik-( rint(delyik/(L*sin(THETA)) )*L*sin(THETA) );

				Rik = sqrt( pow(delxik,2)+pow(delyik,2) );

//				printf("Rik=%f for sites %d and %d\n",Rik,i,connect[i][MOD(j+3,6)]);
//				printf("Bond=%d\n",(j+3)%6);
//				&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		

				CROSSjik=delxij*delyik-delyij*delxik;

				THETAjik=asin(sqrt(pow(CROSSjik,2))/(Rij*Rik));    // david's original small angle approximation
//				THETAjik=sqrt(pow(CROSSjik,2))/(Rij*Rik);		// we will use the sin in energy calculation

				THETAjik_half=THETAjik/2.0;    // theta / 2
				sin_THETAjik_half = sin(THETAjik_half);   // sin (theta/2)

				if (func_flag == 1){
				i_new = i;               // used to write new line

//				if(j == 2){   // this condition only works when allbonds are present
				if(i_prev != i_new){   // this condition is for pbond < 1
					fprintf(fout,"\n "); 		// after three angles (180 degree each initially at each node) new line for new node position
				}

//				fprintf(fout,"%d ", i);  // write bending angles to output file				
//				fprintf(fout,"%d ", j);  // write bending angles to output file				
//				fprintf(fout,"%f ", asin(THETAjik));  // write bending angles to output file
				fprintf(fout,"%f ", asin(CROSSjik/(Rij*Rik)));  // write bending angles to output file
				i_prev = i;             // used to write new line

				}


//				printf("THETAjik=%f\n",THETAjik*180/M_PI);
//				printf("\n");

//				FILB[3*i+j]=kappa2/2*pow(THETAjik,2);
				FILB[3*i+j]=2*kappa2*pow(sin_THETAjik_half,2);   // 2*kappa*sin squared theta/2
		
				Filbend=Filbend+FILB[3*i+j];
//				printf("\n");


			}
			if( (bone==0 || btwo==0) ){
				if (func_flag == 1){
				i_new = i;               // used to write new line

//				if(j == 2){   // this condition only works when allbonds are present
				if(i_prev != i_new){   // this condition is for pbond < 1
					fprintf(fout,"\n "); 		// after three angles (180 degree each initially at each node) new line for new node position
				}

//				fprintf(fout,"%d ", i);  // write bending angles to output file				
//				fprintf(fout,"%d ", j);  // write bending angles to output file				
				fprintf(fout,"%s", "NA ");  // write bending angles to output file
				i_prev = i;             // used to write new line

				}

			}
		}
	}
}

//	ETOT = Eelastic/(UNITC/4)+Ebend/(UNITC/4)+Filbend/(UNITC/4);
//	ETOT = Eelastic/(UNITC/4)+Ebend/(UNITC/4)+Filbend/(UNITC/4)+Eforce/(UNITC/4);

	ETOT = Eelastic/(UNITC/4)+Ecomp/(UNITC/4)+Ebend/(UNITC/4)+Filbend/(UNITC/4);

//	printf("\n %f %f %f %f \n", ETOT, Eelastic/(UNITC/4), Ebend/(UNITC/4), Filbend/(UNITC/4));

	if (func_flag == 1){
//		printf("\n %f %f %f %f \n", ETOT, Eelastic/(UNITC/4), Ebend/(UNITC/4), Filbend/(UNITC/4));
		//old energy ouptup file
//		output_energy_sub(ETOT, Eelastic/(UNITC/4), Filbend/(UNITC/4), Eforce/(UNITC/4), Etension/(UNITC/4), Etension_dip/(UNITC/4), Ecompression/(UNITC/4), Ebuckle/(UNITC/4));   // write energies to output file
		output_energy_sub(ETOT, Eelastic/(UNITC/4), Ecomp/(UNITC/4), Filbend/(UNITC/4), Ebend/(UNITC/4));   // write energies to output file   // Ebend is for arp
		fclose(fout);  // close the output file
		free(filename);

		// writing energy at each node to output file		
		output_energy_node_sub(ES, FILB);   // write energies to output file   // Ebend is for arp
	}


	if (func_flag == 1){
		//strain output file
		output_strain_sub();   // write energies to output file
	}

	return ETOT;
}


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&& DERIVATIVE OF FUNC &&&&&&&&&&&&&&&&&&&&&&&&&&

//This subroutine takes the derivative of func for each free parameter in && i>=L && i<N-L
//the lattice. It returns a dvector(array) of length 2*N with each entry new_topology_conj_grad_618_sans_periodic_Anal-dfunc_SUBLAT.c
//being the derivative of each lattice sites x and y coordinate. 

//********* THIS IS THE ANALYTIC DFUNC **********/

void dfunc( double p[], double xi[]){

	double dx[N],dy[N];
	double deltaX,deltaY,rij,rik,rjk;
	double delxij,delyij,delxik,delyik,delxjk,delyjk;
	double delxmj,delymj,delxlk,delylk;
	double THETAjik,THETAijk,THETAjki,fjik,fijk,fjki,DfdX,DfdY;
	double THETAmji,fmji,CROSSmji,CROSSjik,rmj;
	double THETAlki,flki,CROSSlki,rlk;

	double FSx[N*6],FSy[N*6],FBx[N*6],FBy[N*6];//;,FStot[N][6];
	// double FILAx[N*3],FILAy[N*3];

	int i,j,k,m,bone,btwo,bthree,bfour,arp1,arp2,arp3;
	int tempcon1,tempcon2;	

	int dip_flag = 0;

	int fsign_val, jj;

	int dip_true, mu_flag;  // mu_flag is a flag used to distiguish between compression and buckling

	int row, central_dip;
	int dip_cent_flag,kk, dip_conn_flag;

//printf("Dfunc Indeed!\n");

// keeping ALL boundaries fixed
	for(i=0;i<N;i++){

		dip_true = 0;    // flag for checking if we are at a dipole node

		for (jj=0;jj<(num_dip_unique);jj++){
			if ( (i == dip_node_unique[jj]) ) {
				dip_true = 1;
			}
		}
		if (outer_node[i] == 1){      // outer boundary
			dx[i]=pr[i];
			dy[i]=pr[i+N];
		}
		else{
			dx[i]=p[i];
			dy[i]=p[i+N];
			
		}
	}



for(i=0;i<N;i++){
	xi[i+1]=xi[i+N+1]=0.0;
	for(j=0;j<6;j++){
		
		FSx[i*6+j]=FSy[i*6+j]=FBx[i*6+j]=FBy[i*6+j]=0.0;
		if(j/3==0){
			FILAx[3*i+j]=FILAy[3*i+j]=0.0;
		}
	}

}
//!!!!!!!!!!!!!!!!!COMPUTE DERIVATIVE!!!!!!!!!!!!!!!
					//*********************BOND SPRINGS*****************
if(mu > 0.0){				//**************************************************
	for(i=0;i<N;i++){
		for(j=0;j<6;j++){


			if(occupation[i*6+j]==1){

				deltaX=-(dx[i]-dx[ connect[i][j] ]);
//				deltaX=deltaX-( rint(deltaX/L)*L );

				deltaY=-(dy[i]-dy[ connect[i][j] ]);
//				deltaY=deltaY-( rint(deltaY/(L*sin(THETA)) )*L*sin(THETA) );

//				printf("for site %d & %d deltaX=%f and deltaY=%f\n",i,connect[i][j],deltaX,deltaY);

				rij=sqrt( pow(deltaX,2)+pow(deltaY,2) );

//				delta_bond_len = rij-1.0;
				delta_bond_len = rij - rlen_arr[i*6+j];

//				FSx[6*i+j] = -mu*(rij-1.0)*(deltaX/rij);
//				FSy[6*i+j] = -mu*(rij-1.0)*(deltaY/rij);

				FSx[6*i+j] = -mu*(delta_bond_len)*(deltaX/rij);
				FSy[6*i+j] = -mu*(delta_bond_len)*(deltaY/rij);

//				FSx[6*i+j] = -mu*(rij-1.0)*(deltaX/rij);
//				FSy[6*i+j] = -mu*(rij-1.0)*(deltaY/rij);

//				FSx[6*i+j] = -mu*(rij-1.0)*(deltaX/rij);
//				FSy[6*i+j] = -mu*(rij-1.0)*(deltaY/rij);

//				FStot[i][j]=mu*(rij-1.0);
				
			}
		}
	}
}
				//*********************ARP SPRINGS*****************
				//*************************************************


//*********************The angle ordering for computing the energy of the angle springs
//*********************************************************************************************
/*						          j
							 /\
					    (bthree)rjk /  \ rik (bone)
		       				       /    \
		      				      /      \ 
		      				     k--------i                             
						      rik (btwo)
*/

fjik=0.0;
DfdX=0.0;
DfdY=0.0;
THETAjik=0.0;
delxij=delyij=delxik=delyik=0.0;

				
if(kappa1 > 0.0){
	for(i=0;i<N;i++){
		for(j=0;j<6;j++){
			FBx[6*i+j]=FBy[6*i+j]=0.0;

			bone=occupation[i*6+j]; 
			btwo=occupation[i*6+(j+1)%6];
			bthree=occupation[connect[i][j]*6+(j+2)%6];

//		printf("%d %d %d\n", i,connect[i][j],connect[i][(j+1)%6]);

			arp1=arpoccupation[i*6+j];
			arp2=arpoccupation[connect[i][j]*6+(j+2)%6];
			arp3=arpoccupation[connect[i][(j+1)%6]*6+(j+4)%6];			

			tempcon1=connect[i][j];
			tempcon2=connect[i][(j+1)%6];
//if(i>=L && i<N-L){
			// original commented
// 			if(bone==1 && btwo==1 && arp1==1 ){//Triplet jik //&& j>0 && j<3)//ANISO-exclude all spings except the one between 1 and 2)	  
 			if(bone==1 && btwo==1 ){	  
				DfdX=DfdY=0.0;	
				delxij=delyij=0.0;				
				delxik=delyik=0.0;
				rij=rik=0.0;	
//				&&&&&&&&&&&&&&&& 1st Bond &&&&&&&&&&&&&&&&&			  	
//				&&&&&&&&&&&&&&&&&& BOND ij &&&&&&&&&&&&&&&&&

				delxij = (dx[i]-(dx[ tempcon1 ]));	
//				delxij = delxij-( rint(delxij/L)*L );

				delyij = (dy[i]-(dy[ tempcon1 ]));
//				delyij = delyij-( rint(delyij/(L*sin(THETA)) )*L*sin(THETA) );

				rij = sqrt( pow(delxij,2)+pow(delyij,2) );
//				printf("rij=%f for sites %d and %d\n",rij,i,connect[i][j]);
				
	

//				&&&&&&&&&&&&&&&& 2nd Bond &&&&&&&&&&&&&&&&&&
//				&&&&&&&&&&&&&&&&&& BOND ik &&&&&&&&&&&&&&&&&
	
				delxik = (dx[i]-( dx[ tempcon2 ]));
//				delxik = delxik-( rint(delxik/L)*L );

				delyik = (dy[i]-( dy[ tempcon2 ] ));
//				delyik = delyik-( rint(delyik/(L*sin(THETA)) )*L*sin(THETA) );

				rik = sqrt( pow(delxik,2)+pow(delyik,2) );
//				printf("rik=%f for sites %d and %d\n",rik,i,connect[i][(j+1)%6]);
//				printf("bone=%d and btwo=%d\n",bone,btwo);
//				printf("\n");
//				&&&&&&&&&&&&&& DOT PRODUCT &&&&&&&&&&&&&&&&&				
				fjik=(delxij*delxik+delyij*delyik)/(rij*rik);
				THETAjik=acos(fjik);
			
				DfdX=(delxij+delxik)/(rij*rik)-delxij*fjik/pow(rij,2)-delxik*fjik/pow(rik,2);
				DfdY=(delyij+delyik)/(rij*rik)-delyij*fjik/pow(rij,2)-delyik*fjik/pow(rik,2);
//				&&&&&&&&&&&&&& DERIVATIVE &&&&&&&&&&&&&&&&&
		
				// original commented below
				//FBx[6*i+j]=FBx[6*i+j]+(-kappa1*( (THETAjik-THETA)/sqrt(1-pow(fjik,2)) )*DfdX);
				//FBy[6*i+j]=FBy[6*i+j]+(-kappa1*( (THETAjik-THETA)/sqrt(1-pow(fjik,2)) )*DfdY);

				if (THETAjik < arp_lim){    // checking if two bonds are about to cross beyond 1 degree separation
					FBx[6*i+j]=FBx[6*i+j]+(-kappa1*( (THETAjik-arp_lim)/sqrt(1-pow(fjik,2)) )*DfdX);
					FBy[6*i+j]=FBy[6*i+j]+(-kappa1*( (THETAjik-arp_lim)/sqrt(1-pow(fjik,2)) )*DfdY);
				}
				else{
					FBx[6*i+j]=FBx[6*i+j] + 0;
					FBy[6*i+j]=FBy[6*i+j] + 0;					
				}

			}

			//if(bone==1 && bthree==1 && arp2==1){ //Triplet ijk
			if(bone==1 && bthree==1){ 

//				&&&&&&&&&&&&&&&& 1st Bond &&&&&&&&&&&&&&&&&&		
//				&&&&&&&&&&&&&&&&&& BOND ij &&&&&&&&&&&&&&&&&
				DfdX=DfdY=0.0;	
				delxij=delyij=0.0;				
				delxjk=delyjk=0.0;
				rij=rjk=0.0;	

				delxij = (dx[i]-(dx[ tempcon1 ]));	
				delxij = delxij-( rint(delxij/L)*L );

				delyij = (dy[i]-(dy[ tempcon1 ]));
				delyij = delyij-( rint(delyij/(L*sin(THETA)) )*L*sin(THETA) );

				rij = sqrt( pow(delxij,2)+pow(delyij,2) );
//				printf("rji=%f for sites %d and %d\n",rij,i,connect[i][j]);
//				&&&&&&&&&&&&&&&& 2nd Bond &&&&&&&&&&&&&&&&&&
//				&&&&&&&&&&&&&&&&&& BOND jk &&&&&&&&&&&&&&&&&

				delxjk = (dx[tempcon2]-dx[tempcon1]);
//				delxjk = delxjk-( rint(delxjk/L)*L );
				
				delyjk = (dy[tempcon2]-dy[tempcon1]);
//				delyjk = delyjk-( rint(delyjk/(L*sin(THETA)) )*L*sin(THETA) );
				
				rjk = sqrt( pow(delxjk,2)+pow(delyjk,2));
//				printf("rjk=%f for sites %d and %d\n",rjk,connect[i][j],connect[i][(j+1)%6]);
//				printf("bone=%d and bthree=%d \n",bone,bthree);
//				printf("\n");
//				&&&&&&&&&&&&&& DOT PRODUCT &&&&&&&&&&&&&&&&&
	
				fijk=(delxij*delxjk+delyij*delyjk)/(rij*rjk);
				THETAijk=acos(fijk);

//				&&&&&&&&&&&&&& DERIVATIVE &&&&&&&&&&&&&&&&&	

				DfdX=delxjk/(rij*rjk)-(fijk*delxij)/pow(rij,2);
				DfdY=delyjk/(rij*rjk)-(fijk*delyij)/pow(rij,2);

				// original commented below
				//FBx[6*i+j]=FBx[6*i+j]+(-kappa1*( (THETAijk-THETA)/sqrt(1-pow(fijk,2)) )*DfdX);
				//FBy[6*i+j]=FBy[6*i+j]+(-kappa1*( (THETAijk-THETA)/sqrt(1-pow(fijk,2)) )*DfdY);

				if (THETAijk < arp_lim){    // checking if two bonds are about to cross beyond 1 degree separation
					FBx[6*i+j]=FBx[6*i+j]+(-kappa1*( (THETAijk-arp_lim)/sqrt(1-pow(fijk,2)) )*DfdX);
					FBy[6*i+j]=FBy[6*i+j]+(-kappa1*( (THETAijk-arp_lim)/sqrt(1-pow(fijk,2)) )*DfdY);
				}
				else{
					FBx[6*i+j]=FBx[6*i+j] + 0.0;
					FBy[6*i+j]=FBy[6*i+j] + 0.0;
				}

			}

			//if(btwo==1 && bthree==1 && arp3==1){ // Triplet jki
			if(btwo==1 && bthree==1){ // Triplet jki
				DfdX=DfdY=0.0;
				delxik=delyik=0.0;
				delxjk=delyjk=0.0;
				rik=rjk=0.0;

//				&&&&&&&&&&&&&&&& 1st Bond &&&&&&&&&&&&&&&&&&		
//				&&&&&&&&&&&&&&&&&& BOND ik &&&&&&&&&&&&&&&&&
				delxik = (dx[i]-( dx[ tempcon2 ]));
//				delxik = delxik-( rint(delxik/L)*L );

				delyik = (dy[i]-( dy[ tempcon2 ] ));
//				delyik = delyik-( rint(delyik/(L*sin(THETA)) )*L*sin(THETA) );

				rik = sqrt( pow(delxik,2)+pow(delyik,2) );
//				printf("rik=%f for sites %d and %d\n",rik,i,connect[i][(j+1)%6]);
//				&&&&&&&&&&&&&&&& 2nd Bond &&&&&&&&&&&&&&&&&&
//				&&&&&&&&&&&&&&&&&& BOND jk &&&&&&&&&&&&&&&&&

				delxjk = (dx[tempcon1]-dx[tempcon2]);
//				delxjk = delxjk-( rint(delxjk/L)*L );
				
				delyjk = (dy[tempcon1]-dy[tempcon2]);
//				delyjk = delyjk-( rint(delyjk/(L*sin(THETA)) )*L*sin(THETA) );
				
				rjk = sqrt( pow(delxjk,2)+pow(delyjk,2));
//				printf("rjk=%f for sites %d and %d\n",rjk,connect[i][(j+1)%6],connect[i][j]);
//				printf("btwo=%d and bthree=%d \n",btwo,bthree);
//				printf("\n");
//				&&&&&&&&&&&&&& DOT PRODUCT &&&&&&&&&&&&&&&&&
	
				fjki=(delxik*delxjk+delyik*delyjk)/(rik*rjk);
				THETAjki=acos(fjki);

//				&&&&&&&&&&&&&& DERIVATIVE &&&&&&&&&&&&&&&&&	

				DfdX=delxjk/(rik*rjk)-(fjki*delxik)/pow(rik,2);
				DfdY=delyjk/(rik*rjk)-(fjki*delyik)/pow(rik,2);

				// original commented below
				//FBx[6*i+j]=FBx[6*i+j]+(-kappa1*( (THETAjki-THETA)/sqrt(1-pow(fjki,2)) )*DfdX);
				//FBy[6*i+j]=FBy[6*i+j]+(-kappa1*( (THETAjki-THETA)/sqrt(1-pow(fjki,2)) )*DfdY);

				if (THETAjki < arp_lim){    // checking if two bonds are about to cross beyond 1 degree separation
					FBx[6*i+j]=FBx[6*i+j]+(-kappa1*( (THETAjki-arp_lim)/sqrt(1-pow(fjki,2)) )*DfdX);
					FBy[6*i+j]=FBy[6*i+j]+(-kappa1*( (THETAjki-arp_lim)/sqrt(1-pow(fjki,2)) )*DfdY);
				}
				else{
					FBx[6*i+j]=FBx[6*i+j] + 0.0;
					FBy[6*i+j]=FBy[6*i+j] + 0.0;
				}


			}
//}	
		}
	}
}



			//*********************FILAMENT SPRINGS*************
			//**************************************************
if(kappa2 > 0.0){
	for( i=0;i<N;i++ ){
		for( j=0;j<3;j++ ){
			
			bone=occupation[i*6+j]; 
			btwo=occupation[i*6+(j+3)%6];  
			bthree=occupation[connect[i][j]*6+j];
			bfour=occupation[connect[i][(j+3)%6]*6+(j+3)%6];
	
			fjik=0.0;
			DfdX=0.0;
			DfdY=0.0;
			THETAjik=0.0;
			delxij=delyij=delxik=delyik=0.0;
			fmji=0.0;
			DfdX=0.0;
			DfdY=0.0;
			THETAmji=0.0;
			delxij=delyij=delxmj=delymj=0.0;
			rij=rmj=0.0;
			CROSSmji=0.0;
			CROSSjik=0.0;
			flki=0.0;
			DfdX=0.0;
			DfdY=0.0;
			THETAlki=0.0;
			delxik=delyik=delxlk=delylk=0.0;
			rik=rlk=0.0;
			CROSSlki=0.0;

//			if(  (i<L && connect[i][j]>=N-L)  || (i>=N-L &&  connect[i][(j+3)%6] < L) ){
//				FILAx[3*i+j]=FILAy[3*i+j]=0.0;
//			}else 
			if(bone==1 && btwo==1){  
//				&&&&&&&&&&&&&&&& 1st Bond &&&&&&&&&&&&&&&&&				  	

				delxij = -(dx[i]-(dx[ connect[i][j] ]));	
//				delxij = delxij-( rint(delxij/L)*L );

				delyij = -(dy[i]-(dy[ connect[i][j] ]));
//				delyij = delyij-( rint(delyij/(L*sin(THETA)) )*L*sin(THETA) );

				rij = sqrt( pow(delxij,2)+pow(delyij,2) );

//				&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		
//				&&&&&&&&&&&&&&&& 2nd Bond &&&&&&&&&&&&&&&&&

				delxik = -(dx[i]-( dx[ connect[i][(j+3)%6] ] ));
//				delxik = delxik-( rint(delxik/L)*L );

				delyik = -(dy[i]-( dy[ connect[i][(j+3)%6] ] ));
//				delyik = delyik-( rint(delyik/(L*sin(THETA)) )*L*sin(THETA) );

				rik = sqrt( pow(delxik,2)+pow(delyik,2) );

//				&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		


				CROSSjik=-(delxij*delyik-delyij*delxik);
				fjik=sqrt(pow(CROSSjik,2))/(rij*rik);
				THETAjik=asin(fjik);
				if(fabs(CROSSjik)>0.0){
					DfdX=(delyik-delyij)*CROSSjik/(fabs(CROSSjik)*rij*rik)-fjik*delxij/pow(rij,2)-fjik*delxik/pow(rik,2);
				
					DfdY=(delxij-delxik)*CROSSjik/(fabs(CROSSjik)*rij*rik)-fjik*delyij/pow(rij,2)-fjik*delyik/pow(rik,2);

//					FILAx[3*i+j]=FILAx[3*i+j]+kappa2*(THETAjik)/( sqrt(1-pow(fjik,2)) )*DfdX;
//					FILAy[3*i+j]=FILAy[3*i+j]+kappa2*(THETAjik)/( sqrt(1-pow(fjik,2)) )*DfdY;

// trying the derivative of sin theta squared and not theta squared:
					FILAx[3*i+j]=FILAx[3*i+j]+kappa2*((fjik)/( sqrt(1-pow(fjik,2)) ))*DfdX;
					FILAy[3*i+j]=FILAy[3*i+j]+kappa2*((fjik)/( sqrt(1-pow(fjik,2)) ))*DfdY;
				}
			}
//			}else 
			if(bone==1 && bthree==1){// && btwo != 1){ 
				DfdX=0.0;
				DfdY=0.0;
				delxij=delyij=delxmj=delymj=0.0;
				rij=rmj=0.0;

				tempcon1=connect[i][j];
				tempcon2=connect[connect[i][j]][j];
//				&&&&&&&&&&&&&&&& 1st Bond ij &&&&&&&&&&&&&&&&&				  	

				delxij = -(dx[i]-(dx[ tempcon1 ]));	
//				delxij = delxij-( rint(delxij/L)*L );

				delyij = -(dy[i]-(dy[ tempcon1 ]));
//				delyij = delyij-( rint(delyij/(L*sin(THETA)) )*L*sin(THETA) );

				rij = sqrt( pow(delxij,2)+pow(delyij,2) );

//				&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		
//				&&&&&&&&&&&&&&&& 2nd Bond jk &&&&&&&&&&&&&&&&&

				delxmj = (dx[tempcon1]-( dx[ tempcon2 ] ));
//				delxmj = delxmj-( rint(delxmj/L)*L );

				delymj = (dy[tempcon1]-( dy[ tempcon2 ] ));
//				delymj = delymj-( rint(delymj/(L*sin(THETA)) )*L*sin(THETA) );

				rmj = sqrt( pow(delxmj,2)+pow(delymj,2) );

//				&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		

				CROSSmji=delxij*delymj-delyij*delxmj;
				fmji=sqrt(pow(CROSSmji,2))/(rij*rmj);
				THETAmji=asin(fmji);

				if(fabs(CROSSmji)>0.0){			
			
					DfdX=-delymj*CROSSmji/(fabs(CROSSmji)*rij*rmj)-delxij*fmji/(pow(rij,2));
					DfdY=delxmj*CROSSmji/(fabs(CROSSmji)*rij*rmj)-delyij*fmji/(pow(rij,2));

//					FILAx[3*i+j]=FILAx[3*i+j]+kappa2*(THETAmji)/( sqrt(1-pow(fmji,2)) )*DfdX;
//					FILAy[3*i+j]=FILAy[3*i+j]+kappa2*(THETAmji)/( sqrt(1-pow(fmji,2)) )*DfdY;

// trying the derivative of sin theta squared and not theta squared:
					FILAx[3*i+j]=FILAx[3*i+j]+kappa2*((fmji)/( sqrt(1-pow(fmji,2)) ))*DfdX;
					FILAy[3*i+j]=FILAy[3*i+j]+kappa2*((fmji)/( sqrt(1-pow(fmji,2)) ))*DfdY;
//					printf("dx=%f dy=%f\n",FILAx[3*i+j],FILAy[3*i+j]);
				}
			}
//			}else 
			if(btwo==1 && bfour==1){// && bone !=1){  	 

				DfdX=0.0;
				DfdY=0.0;
				delxik=delyik=delxlk=delylk=0.0;
				rik=rlk=0.0;

				tempcon1=connect[i][(j+3)%6];
				tempcon2=connect[connect[i][(j+3)%6]][(j+3)%6];
//				&&&&&&&&&&&&&&&& 1st Bond ij &&&&&&&&&&&&&&&&&				  	

				delxik = -(dx[i]-(dx[ tempcon1 ]));	
//				delxik = delxik-( rint(delxik/L)*L );

				delyik = -(dy[i]-(dy[ tempcon1 ]));
//				delyik = delyik-( rint(delyik/(L*sin(THETA)) )*L*sin(THETA) );

				rik = sqrt( pow(delxik,2)+pow(delyik,2) );

//				&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		
//				&&&&&&&&&&&&&&&& 2nd Bond jk &&&&&&&&&&&&&&&&&

				delxlk = (dx[tempcon1]-( dx[ tempcon2 ] ));
//				delxlk = delxlk-( rint(delxlk/L)*L );

				delylk = (dy[tempcon1]-( dy[ tempcon2 ] ));
//				delylk = delylk-( rint(delylk/(L*sin(THETA)) )*L*sin(THETA) );

				rlk = sqrt( pow(delxlk,2)+pow(delylk,2) );

//				&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&		


				CROSSlki=-(delyik*delxlk-delxik*delylk);
				flki=sqrt(pow(CROSSlki,2))/(rik*rlk);
				THETAlki=asin(flki);

				if(fabs(CROSSlki) > 0.0){
					DfdX=-delylk*CROSSlki/(fabs(CROSSlki)*rik*rlk)-delxik*flki/(pow(rik,2));
					DfdY=delxlk*CROSSlki/(fabs(CROSSlki)*rik*rlk)-delyik*flki/(pow(rik,2));

//					FILAx[3*i+j]=FILAx[3*i+j]+kappa2*(THETAlki)/( sqrt(1-pow(flki,2)) )*DfdX;
//					FILAy[3*i+j]=FILAy[3*i+j]+kappa2*(THETAlki)/( sqrt(1-pow(flki,2)) )*DfdY;

// trying the derivative of sin theta squared and not theta squared:
					FILAx[3*i+j]=FILAx[3*i+j]+kappa2*((flki)/( sqrt(1-pow(flki,2)) ))*DfdX;
					FILAy[3*i+j]=FILAy[3*i+j]+kappa2*((flki)/( sqrt(1-pow(flki,2)) ))*DfdY;
					//printf("dx=%f dy=%f\n",FILAx[3*i+j],FILAy[3*i+j]);
				}
			}
		}
	}
}


	for(i=0;i<N;i++){

		xi[i+1]=0.0;
		xi[i+N+1]=0.0;

		// xi above stores forces but forces on boundary need to be zero
		// so I created fxi which is the same as xi but stores real forces on boundary nodes
		// fxi is not used in CG method. It is just for analysis of boundary forces
		fxi[i]=0.0;
		fxi[i+N]=0.0;

		// forces only work on rows between the top and bottom rows like in david's original code: this is now for nodes that are not on 
		// any boundary. ALL boundary nodes are excluded!
//		if(i>=L && i<N-L && i%L != 0 && i%L != L-1 && dip_true == 0){		
//		if(i>=L && i<N-L && i%L != 0 && i%L != L-1){		// rectangular boundary fixed 
		if (outer_node[i] == 0){      // outer circular boundary
			for(j=0;j<6;j++){

				xi[i+1]=xi[i+1]+FSx[6*i+j]+FBx[6*i+j];
				xi[i+N+1]=xi[i+N+1]+FSy[6*i+j]+FBy[6*i+j];

				if(j/3==0){
						
					xi[i+1]=xi[i+1]+FILAx[3*i+j];
					xi[i+N+1]=xi[i+N+1]+FILAy[3*i+j];
					//printf("dx=%f dy=%f\n",xi[i+1],xi[i+N+1]);
				}
			}
		}

			for(j=0;j<6;j++){

				fxi[i]=fxi[i]+FSx[6*i+j]+FBx[6*i+j];
				fxi[i+N]=fxi[i+N]+FSy[6*i+j]+FBy[6*i+j];

				if(j/3==0){
						
					fxi[i]=fxi[i]+FILAx[3*i+j];
					fxi[i+N]=fxi[i+N]+FILAy[3*i+j];
					//printf("dx=%f dy=%f\n",xi[i+1],xi[i+N+1]);
				}
			}
		}



}




//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&& MAIN &&&&&&&&&&&&&&&&&&&&&&&&&&&&

//int main(/* int argc, char **argv*/ )
int main( int argc, char **argv )
{ 

	int i,j,xpos,ypos,r,k,ms, row, col, column, jj; 
	int iter,t,MAXREALi,pstep,tstep;
	double prob,pfila,parp,mumu;
	double ftol,fret,ass,delta,gmmaXY,gmmaXX,smod,kapmod,sigmod;
	double LdivLam,perc,occusum,percrat,freqbond,freqarp,arpsum,realP,seed,factor,GAMMA1,GAMMA2;
	double angle; // this is the angle between x axis and two nodes of a dipole for the single dipole case (used to create dipoles that are not on x axis)

	int d_i;
	int sign_val = 1;  // to alternate sign of addition or subtraction of perturbation of points
	int dip_flag = 0;
	int old_times; // used only if hysterisis is true. Simulation will run once towards max perturbation and then back to zero
	double pert_val, dx_val;
	int fsign_val; 		// sign of force that alternates between positive and negative
	int kk;
	int num_node_dip;

	double move_119x, move_119y, move_120x, move_120y, move_134x, move_134y, move_135x, move_135y;
	double move_136x, move_136y, move_151x, move_151y, move_152x, move_152y;

	int dip_true;

	int loc_flag;   // for checking that the dipole is not at the edge

	double dx_node_center, dy_node_center;   // to find whether we are in the inner or outer circle

	float rand_val;

	double dx_center, dy_center, dist_center;
	double dr_center = 0.5;

	iter=0;	
	ftol=TOL;

	gmmaXY=0.0;
	gmmaXX=0.0;
	smod=0.0;
	prob=0.0;

	pstep=0;
	MAXREALi=1;
	times=0;
	
	double C44[pstep+1],ddmod[times],sigma[times],F[2*N+1];//test[2*N];
	

//------------------------------------------------------------------------------------	
//                             BASH CODE below
//------------------------------------------------------------------------------------	
	printf("number of arguments: %d\n", argc);
    for (i=0; i < argc; i++)
    {
    	printf("argument %d: (%s)\n", i, argv[i]);
    }


   	if (argc < 2)
	{
		printf("Need input file.\n");
		return(0);
	}

	char paramfname[255], string [255];
	double mu_test, pbond_test;
//	strcpy(paramfname, argv[1]);
//	printf("\n paramter file name: %s \n", paramfname);
	strcpy(arg_txt,argv[1]);
	sscanf(argv[1], "%ld", &ranseed1);
	printf("ranseed1_bash: %ld \n", ranseed1);
	sscanf(argv[2], "%ld", &ranseed2);
	printf("ranseed2_bash: %ld \n", ranseed2);
	sscanf(argv[3], "%d", &diff_network1);
	printf("diff_network1: %d \n", diff_network1);
	sscanf(argv[4], "%d", &diff_network2);
	printf("diff_network2: %d \n", diff_network2);
	sscanf(argv[5], "%lf", &kappa2);
	printf("kappa_bash: %f \n", kappa2);
	sscanf(argv[6], "%lf", &mu_c);
	printf("mu_c_bash: %f \n", mu_c);
	sscanf(argv[7], "%lf", &PBOND);
	printf("pbond_bash: %f \n", PBOND);
	sscanf(argv[8], "%d", &num_center);
	printf("num_center_bash: %d \n", num_center);
	sscanf(argv[9], "%d", &srand_arg);
	printf("srand val: %d \n", srand_arg);

//------------------------------------------------------------------------------------	

	
//	num_center = 1;

	num_dip = 6*num_center;
	num_node = 2*num_dip;
	num_node_dip = 12;   // number of nodes per center of contraction: 2*6 dips = 12

	dip_node = calloc(num_node, sizeof(int));
	Ff = calloc(num_node, sizeof(double));
	loc_center = calloc(num_center, sizeof(int));



//************************************************************
//**********************LATTICE POSITIONS (x,y)*************** 

	for (i=0;i<N;i++) {
	
//!!!!!!!!!!!!!!!!On the square lattice!!!!!!!!!!!!!!!!!!!!!!!
		xpos=i%L;	
		ypos=i/L;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!On the triangular lattice!!!!!!!!!!!!!!!!!!!

//We will create a lattice that is L*sin(THETA) wide in a box that is L wide
//In order to keep the bod length constrained to 1 we need to "squash" the 
//lattice in the y-coordinate


		y[i]=(double)(ypos)*sin(THETA); 
		x[i]=xpos;
								//	<--------------L------------> ^
//Shift the odd labeled rows. if (ypos=odd)			//	@	@	@ 	@     |
								//	    @	    @       @  	    @ |	
		if( ypos%2 == 1 ){				//	@	@	@ 	@     L*sin(THETA)
								//	    @	    @       @       @ |
			x[i]=x[i]-cos(THETA); 			//	@	@	@ 	@     |
		}						//				      V
	
		if (i == 4053){
			printf("4053 x: %0.10f y: %0.10f ypos: %d \n",x[4053], y[4053], ypos);
		}
		if (i == 4095){
			printf("4095 x: %0.10f y: %0.10f ypos: %d \n",x[4095], y[4095], ypos);
		}
	}							//					
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//************************************************************
//*******************TOPOLOGY AND PERIODICITY*****************


// 	DIAGRAM 1			DIAGRAM 2
//	(even)				(odd)

/*      (2)     (1)		       (1)     (2)
          \      |			|      /
           \     |			|     /
            \    |			|    /
             \   |			|   /
              \  |			|  /	
               \ | 			| /
   (0)___________[ms]   (0)____________[ms] */   

j=0;
	for(j=0;j<L;j++){

		for(i=0;i<L;i++){

//!!!!!!!!!!!!!!!!ms is the site index!!!!!!!!!!!!!!!!!!!!!!!!! 
//!!!!!!(i) is xpos and (j) is ypos on the square lattice!!!!!!

			ms=i+j*L;

//!!!!!This is for the periodic BC's in the x-direction(BOTH DIAGRAMS)!!!!!

			connect[ms][0] = MOD(i-1,L)+j*L; 
			connect[ms][3] = MOD(i+1,L)+j*L;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!! DIAGRAM 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
			if(j%2 == 0){

//			
//Periodic BC's for the y-direction,j is the bond number	
				connect[ms][1] = i+MOD(j-1,L)*L;	//For the ypos connection was 2
				connect[ms][5] = i+MOD(j+1,L)*L;	//The reflection was 4			

//If ypos is even the diagonal connection is up and to the left
				connect[ms][2] = MOD(i+1,L)+MOD(j-1,L)*L;	//For the Diagonal connection was 1
//				connect[ms][5] = MOD(i-1,L)+MOD(j+1,L)*L;	//Its relection was 5
				connect[ms][4] = MOD(i+1,L)+MOD(j+1,L)*L;
			}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
//!!!!!!!!!!!!!!!!!!!!!!! DIAGRAM 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
			if(j%2 == 1){ 
//Periodic BC's for the y-direction,j is the bond number				
				connect[ms][2] = i+MOD(j-1,L)*L; //was 1
				connect[ms][4] = i+MOD(j+1,L)*L; //was 5
//If ypos is odd the diagonal connection is up and to the right
				connect[ms][1] = MOD(i-1,L)+MOD(j-1,L)*L; //was i+1 and j-1
				connect[ms][5] = MOD(i-1,L)+MOD(j+1,L)*L;
			}
		}	
	}

// after the neighbor list is populated, we can start by setting up dipoles for isotropic contraction.
//	loc_center = 2016;

	int boundary_sep = 20;    // least separation from the boundary
	int dip_sep;
	int r_in = 12;

	if (num_center > 50){
		dip_sep = 2.5;       // least separation from other dipoles:: for trying 100 dipoles
	}
	else{
		dip_sep = 2.5;       // least separation from other dipoles:: original for all data untill March 19 2024
	}

	int dx_dip, dy_dip;

	// first define center of the circular network
//	int center = 2016;
	int center = 2080;
	int r_out = 25;

	if (L == 128){
		center = 8383;
		r_in = 24;
		r_out = 50;
	}
	if (L == 65){
		center = 2113;
	}

	double dx_dip_center, dy_dip_center;

//	srand(111); 

	srand(srand_arg); 

	for (i = 0; i < num_center; i++){
		loc_flag = 0;

		while(loc_flag == 0){    // checking if the dipole is not at the edge
		loc_center[i] = rand()%(((N-1)+1)-0);   // N-1 is max, 0 is min; +1 needed
		if (num_center == 1 && L == 16){
//		if (num_center == 1){
			loc_center[i] = 135;   // for one contraciton point, place it in the center of the box			
		}
		else if (num_center == 1 && L == 64){
//			loc_center[i] = 2016;   // for one contraciton point, place it in the center of the box			
//			loc_center[i] = 2080;   // for one contraciton point, place it in the center of the box			
		}
		else if (num_center == 1 && L == 65){
			loc_center[i] = 2113;   // for one contraciton point, place it in the center of the box	
//			center = loc_center[i];	
		}


//		row = (int)loc_center[i]/L;
//		column = loc_center[i]%L; 
		dx_dip_center = x[center] - x[loc_center[i]]; 
		dy_dip_center = y[center] - y[loc_center[i]]; 

		if (sqrt((dx_dip_center*dx_dip_center + dy_dip_center*dy_dip_center)) > r_in) loc_flag = 0;

		// dipoles cannot be near boundaries
//		if ((row <= boundary_sep) || (row >= (L-boundary_sep-1)) || (column <= boundary_sep) || (column >= (L-boundary_sep-1))) loc_flag = 0;   // no dipoles in top or bottom 10 rows
//		else if ((column <= boundary_sep) || (column >= (L-boundary_sep-1))) loc_flag = 0;   // no dipoles near the boundaries
		else loc_flag = 1;

		// start from second dipole center and compare with every center defined before
//		if (i >= 1){
		if ((i >= 1) && (loc_flag != 0)){
			for (k=0;k<i;k++){

				dx_dip = x[loc_center[i]] - x[loc_center[k]]; 
				dy_dip = y[loc_center[i]] - y[loc_center[k]]; 

				if (sqrt((dx_dip*dx_dip + dy_dip*dy_dip)) <= dip_sep ) loc_flag = 0;

				if (loc_flag == 0){
					printf(" not accepting %d because it is too close to %d \n", loc_center[i], loc_center[k]);
				}
			}
		}

		printf("\n location: %d  row: %d   column: %d   loc_flag: %d \n", loc_center[i], row, column, loc_flag);
		}
	}

	int test_flag = 0;   // this flag to test the e_bend and e_stretch by manually placing dipoles in the network.
// checking the code by manually putting dipoles at different distances
	if (num_center == 2 && L == 64 && test_flag == 1){
		printf("\n entering manual liftoff! Location specified! \n");
		loc_center[0] = 2016;   // for one contraciton point, place it in the center of the box			
		loc_center[1] = 2019;   // for one contraciton point, place it in the center of the box			
	}
	if (num_center == 1 && L == 64 && test_flag == 1){
		printf("\n entering manual liftoff! Location specified! \n");
//		loc_center[0] = 966;   // for one contraciton point, place it in the center of the box			
		loc_center[0] = 1987;   // for one contraciton point, place it in the center of the box			
	}



	// creating an array of unique nodes
	int counter_dip = 0;
	num_dip_unique = 7*num_center;
	dip_node_unique = calloc(num_dip_unique, sizeof(int));  // 7 because each isotropic dipole has one node in center 6 nodes around it

	printf("\n Printing nodes of dipoles \n");
	for (j = 0; j< num_center; j++)
	{
		kk=0;

		dip_node_unique[counter_dip] = loc_center[j];
		printf("\n dip is:: %d ",dip_node_unique[counter_dip]);
		counter_dip = counter_dip + 1;

//		printf("\n j = %d  loc_center = %d \n",j,loc_center[j]);
		for (i=0;i<num_node_dip;i++)
		{

			if (i%2 == 0) dip_node[i+j*num_node_dip] = loc_center[j];
//			if (i%2 == 0) dip_node[i] = loc_center;
			else if (i%2 == 1)
			{
				dip_node[i+j*num_node_dip] = connect[loc_center[j]][kk];
//				dip_node[i] = connect[loc_center][kk];
//				printf("\n i: %d kk: %d      %d \n ",i,kk,connect[loc_center[j]][kk]);
				kk = kk+1;

				dip_node_unique[counter_dip] = dip_node[i+j*num_node_dip];     // adding unique nodes to array of dipoles
				printf("\n dip at %d is : %d ",counter_dip, dip_node_unique[counter_dip]);
				counter_dip = counter_dip + 1;
			}
//			printf("\n node[%d] = %d \n",i,dip_node[i+j*num_node_dip]);
		}

	}

	for (i = 0; i < num_dip_unique; i++){
		printf("\n num_dip_unique = %d \n",dip_node_unique[i]);
	}

	printf("\n   -------   \n");

	int temp;
	printf("\n num nodes: %d \n", num_node);
	for (i=0;i<(num_node-1);i=i+2)
	{
//		printf("\n %d  %d \n ",dip_node[i],dip_node[i+1]);
		if (x[dip_node[i+1]] < x[dip_node[i]]){
//			printf("\n %d  %d \n ",dip_node[i],dip_node[i+1]);
			temp = dip_node[i+1];
			dip_node[i+1] = dip_node[i];
			dip_node[i] = temp;
		}
	}

	printf("\n number of dipoles: %d \n", num_dip);
	printf("\n number of nodes: %d \n", num_node);
//	printf("\n node[%d] = %d \n",16,dip_node[16]);
//	printf("\n node[%d] = %d \n",17,dip_node[17]);


for ( i=0; i<N; i++ ) {
	for (j=0;j<(num_node);j++){
		if (i == dip_node[j]){			
			Ff[j] = 0;
//			printf("\n Inital extra force: %f at %d \n", Ff[j], i);
		}
	}
} 	

	for (i=0;i<N;i++) { 
			subpr[i]=pr[i]=x[i];//	Need to reset positions for nex realization
			subpr[i+N]=pr[i+N]=y[i];//
	}


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& PERCOLATION ROUTINE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& COUPLED BOND AND ARP SPRING PERCOLATION &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//ranseed1 = 667720;  // 667720; 167720; 712240; 537206; 167720; 267720; extra networks: 367725; 127020; 507225; 399020
//ranseed2 = 601210;  // 601210; 101210; 203915; 701210; 101210; 201210; extra networks: 927820; 441039; 102270; 800099 

RANSEED(ranseed1,ranseed2);//  556443,238888;999999,4545454;17765674,9087654 new network 2 dip case

//RANSEED(537206,701210);//  556443,238888;999999,4545454;17765674,9087654
//RANSEED(667720,601210);//  556443,238888;999999,4545454;17765674,9087654 new network 2 dip case

//RANSEED(422041,339098);//  556443,238888;999999,4545454;17765674,9087654
//RANSEED(152240,193915);//  556443,238888;999999,4545454;17765674,9087654
//RANSEED(712240,203915);//  556443,238888;999999,4545454;17765674,9087654
//RANSEED(420149,301961);//  556443,238888;999999,4545454;17765674,9087654 new network 2 dip case

//RANSEED(527149,301961);//  for two dipoles on y axis for clusters
//RANSEED(827450,391961);//  for two dipoles on y axis for clusters


prob=PBOND; // 0.64 defined above (comment by abhinav)
//FILE * dq;
//dq=fopen( "L64-P-0.41_NA_kap1-1e-6_042511_TOL-10e-8.dat","a" );

for(k=0;k<=pstep;k++){
	C44[k]=0.0;
	for(r=0;r<MAXREALi;r++){	
		for (i=0;i<N;i++) { 
//			if (i==i_p) printf("1 subpr[%d] = %f \n", i, subpr[i]);
			subpr[i]=pr[i]=x[i];//	Need to reset positions for nex realization
			subpr[i+N]=pr[i+N]=y[i];//
//			if (i==i_p) printf("2 subpr[%d] = %f \n", i, subpr[i]);
		}

//****************************************************************
//**************************RANDOMIZE*****************************

//	RANSEED(time(NULL),time(NULL)%5);

	int temp_count = 0;  // will be used to cound how big the boundary array should be :0
	int inner_count = 0;  // counts the number of inner nodes
	int outer_all_count = 0;  // all nodes outside inner nodes

	for(i=0;i<N;i++){
		dx_center = x[center] - x[i];
		dy_center = y[center] - y[i];
		dist_center = sqrt(dx_center*dx_center + dy_center*dy_center);

		if (dist_center <= r_in){
//			printf("\n Node number %d is inside the inner circle of radius! %d", i, r_in);
			inner_count += 1;
		}
		else{
//			printf("\n Node number %d is outside the inner circle of radius! %d", i, r_in);
			outer_all_count += 1;
		}
	}

	printf("\n We have a total of %d inner nodes", inner_count);
	// using temp count to define the aray size for boundary nodes.
	int inner_nodes[inner_count];
	int outer_all_nodes[outer_all_count];

	inner_count = 0;
	outer_all_count = 0;

	for(i=0;i<N;i++){
		dx_center = x[center] - x[i];
		dy_center = y[center] - y[i];
		dist_center = sqrt(dx_center*dx_center + dy_center*dy_center);

		if (dist_center <= (r_in)){
			inner_nodes[inner_count] = i;
			inner_count += 1;
		}
		else{
			outer_all_nodes[outer_all_count] = i;
			outer_all_count += 1;			
		}
	}

	float prob_eff;
	int num_bonds_p1_eff = L*L*3 - 4*L;
	printf(" \n Effectively for p = 1 network, the total number of bonds are: %d \n", num_bonds_p1_eff);

	int desired_num_bonds = prob*num_bonds_p1_eff;


	if (PBOND < 1.0)
	{
		prob_eff = desired_num_bonds/(L*L*3.0);		
	}
	else{
		prob_eff = PBOND;
	}

	prob_eff = PBOND;

	printf("The desired number of bonds is: %d \n", desired_num_bonds);
	printf("The effective probability is: %f \n", prob_eff);

	for(i=0;i<inner_count;i++){
		//RANSEED(ranseed1,ranseed2);
		for(j=0;j<3;j++){
//			if( (rand_double() < prob) ){					
			if( (rand_double() < prob_eff) ){					
				occupation[inner_nodes[i]*6+j]=1;
				occupation[connect[inner_nodes[i]][j]*6+(j+3)%6]=1;

			}else{
					occupation[inner_nodes[i]*6+j]=0;
					occupation[connect[inner_nodes[i]][j]*6+(j+3)%6]=0;
			}
		}
	}

	RANSEED(diff_network1,diff_network2);   // re-seeding the ranseed
	for(i=0;i<outer_all_count;i++){
		for(j=0;j<3;j++){
//			if( (rand_double() < prob) ){					
			if( (rand_double() < prob_eff) ){					
				occupation[outer_all_nodes[i]*6+j]=1;
				occupation[connect[outer_all_nodes[i]][j]*6+(j+3)%6]=1;

			}else{
					occupation[outer_all_nodes[i]*6+j]=0;
					occupation[connect[outer_all_nodes[i]][j]*6+(j+3)%6]=0;
			}
		}
	}

// making the boundaries - clamped - so deleting some bonds that were required for PBD condition
	for(i=0;i<N;i++){
			row = (int)i/L;
			column = i%L; 
			if ((row % 2 == 0) && (column == 0)){    // removing 0th bond from even rows on left edge so j = 0 below
				occupation[i*6+0] = 0;
				occupation[connect[i][0]*6+(0+3)%6]=0;
			}
			if ((row % 2 != 0) && (column == 0)){
				occupation[i*6+0] = 0;
				occupation[connect[i][0]*6+(0+3)%6]=0;
				occupation[i*6+1] = 0;
				occupation[connect[i][1]*6+(1+3)%6]=0;
			}

			if ((row % 2 == 0) && (column == L-1)){
				occupation[i*6+2] = 0;
				occupation[connect[i][2]*6+(2+3)%6]=0;
			}
			if (row == 0){
				occupation[i*6+1] = 0;
				occupation[connect[i][1]*6+(1+3)%6]=0;
				occupation[i*6+2] = 0;
				occupation[connect[i][2]*6+(2+3)%6]=0;
			}

	}

// adding bonds to the central dipole node - special node
for (i=0;i<num_center;i++){
	for(j=0;j<6;j++){
		occupation[loc_center[i]*6+j]=1;
 		occupation[connect[loc_center[i]][j]*6+(j+3)%6]=1;
 	}
}

//  edited to enforce radial bonds going out of the dipole node in the matrix
 for (i=0;i<N;i++){
 	for(j=0;j<6;j++){
 		for (k=0;k<num_center;k++){

			if ((i == loc_center[k] - 1) && (j == 0)){
				occupation[i*6+j]=1;
		 		occupation[connect[i][j]*6+(j+3)%6]=1;
			}
			else if ((i == loc_center[k] + 1) && (j == 3)){
				occupation[i*6+j]=1;
		 		occupation[connect[i][j]*6+(j+3)%6]=1;
			}

 			row = (int)loc_center[k]/L;
 			if (row%2 == 0){

				if ((i == loc_center[k] - L) && (j == 1)){
					occupation[i*6+j]=1;
			 		occupation[connect[i][j]*6+(j+3)%6]=1;
				}
				else if ((i == loc_center[k] - L + 1) && (j == 2)){
					occupation[i*6+j]=1;
			 		occupation[connect[i][j]*6+(j+3)%6]=1;
				}
				else if ((i == loc_center[k] + L + 1) && (j == 4)){
					occupation[i*6+j]=1;
			 		occupation[connect[i][j]*6+(j+3)%6]=1;
				}
				else if ((i == loc_center[k] + L) && (j == 5)){
					occupation[i*6+j]=1;
			 		occupation[connect[i][j]*6+(j+3)%6]=1;
				}
 			}

			else{
				if ((i == loc_center[k] - L - 1) && (j == 1)){
					occupation[i*6+j]=1;
			 		occupation[connect[i][j]*6+(j+3)%6]=1;
				}
				else if ((i == loc_center[k] - L) && (j == 2)){
					occupation[i*6+j]=1;
			 		occupation[connect[i][j]*6+(j+3)%6]=1;
				}
				else if ((i == loc_center[k] + L) && (j == 4)){
					occupation[i*6+j]=1;
			 		occupation[connect[i][j]*6+(j+3)%6]=1;
				}
				else if ((i == loc_center[k] + L - 1) && (j == 5)){
					occupation[i*6+j]=1;
			 		occupation[connect[i][j]*6+(j+3)%6]=1;
				}
			}	
		}
	}
}

if (PBOND != 1) {
	//  and removing hexagonal bonds for the dipoles
	 for (i=0;i<N;i++){
	 	for(j=0;j<6;j++){
	 		for (k=0;k<num_center;k++){
	 			row = (int)loc_center[k]/L;

	 			if (row%2 == 0){

					if ((i == loc_center[k] - 1) && (j == 2)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - 1) && (j == 4)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + 1) && (j == 1)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + 1) && (j == 5)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}

					else if ((i == loc_center[k] - L) && (j == 5)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - L) && (j == 3)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - L + 1) && (j == 0)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - L + 1) && (j == 4)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + L + 1) && (j == 2)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + L + 1) && (j == 0)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + L) && (j == 3)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + L) && (j == 1)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
	 			}

				else{
					if ((i == loc_center[k] - 1) && (j == 2)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - 1) && (j == 4)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + 1) && (j == 1)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + 1) && (j == 5)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - L - 1) && (j == 5)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - L - 1) && (j == 3)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - L) && (j == 0)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] - L) && (j == 4)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + L) && (j == 2)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + L) && (j == 0)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + L - 1) && (j == 3)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
					else if ((i == loc_center[k] + L - 1) && (j == 1)){
						occupation[i*6+j]=0;
				 		occupation[connect[i][j]*6+(j+3)%6]=0;
					}
				}	
			}
		}
	}
}


// adding bonds to the central dipole node - special node
for (i=0;i<num_center;i++){
	for(j=0;j<6;j++){
//		if (occupation[loc_center[i]*6+j] == 0)
//		{
			occupation[loc_center[i]*6+j]=1;
	 		occupation[connect[loc_center[i]][j]*6+(j+3)%6]=1;
//		}
 	}
}


// finding the coordination number
	for(i=0;i<N;i++){
		for(j=0;j<6;j++){
			if (occupation[i*6+j]==1){
				coord_num[i] = coord_num[i]+1;
			}
		}
	}


int count, n3, n4, n5, nRandonNumber, jjj, jval, nodeval;
int boundary_flag = 0;

//***************************************************************************************************************************
// the following code sets each node to have at least a connectivity of 3 or 2
	srand(1);

//***************************************************************************************************************************
	// removing dangling bonds
	int dangle_return = dangle();
	printf("\n Dangle return is %d \n", dangle_return);


// removing the dangling bonds by creating another bond in the code below.
	int top_left = (int)((L-1)*L);
	int bottom_right = L-1;


//******************************************************************************************************

	printf("\n Now setting up the circular boundary \n");

	temp_count = 0;  // will be used to cound how big the boundary array should be :0
	int outer_ring_count = 0;

	for(i=0;i<N;i++){
		dx_center = x[center] - x[i];
		dy_center = y[center] - y[i];

		dist_center = sqrt(dx_center*dx_center + dy_center*dy_center);

//		printf("\n For node %d the distance from center %d is %f", i, center, dist_center);
		if ((dist_center >= (r_out-dr_center)) && (dist_center < (r_out+dr_center))){
//			printf("\n Node number %d is on the outer circle of radius! %d", i, r_out);
			temp_count += 1;
		}

		if (dist_center >= (r_out-dr_center)){
			outer_node[i]= 1;
		}

		// counitng the nodes in the outer ring where the network will be changed
		if ((dist_center < (r_out-dr_center)) && (dist_center > r_in)){
				outer_ring_count += 1;
		}
	}

	printf("\n We have a total of %d boundary nodes", temp_count);
	// using temp count to define the aray size for boundary nodes.
	int boundary_nodes[temp_count];
	int outer_ring_nodes[outer_ring_count];

	temp_count = 0;   // counts boundary nodes
	outer_ring_count = 0;

	for(i=0;i<N;i++){
		dx_center = x[center] - x[i];
		dy_center = y[center] - y[i];

		dist_center = sqrt(dx_center*dx_center + dy_center*dy_center);

//		printf("\n For node %d the distance from center %d is %f", i, center, dist_center);
		if ((dist_center >= (r_out-dr_center)) && (dist_center < (r_out+dr_center))){

			boundary_nodes[temp_count] = i;
			temp_count += 1;
		}

		// filling an array that stores the nodes lying in the outer ring
		if ((dist_center < (r_out-dr_center)) && (dist_center > r_in)){
			outer_ring_nodes[outer_ring_count] = i;
			outer_ring_count += 1;
		}
	}

	printf("\n The number of nodes in the outer ring are: %d \n", outer_ring_count);

	// connecting the boundary nodes to each other::
	for(i=0;i<temp_count;i++){
		for (j=0;j<6;j++){
			for(k=0;k<temp_count;k++){
				if (connect[boundary_nodes[i]][j] == boundary_nodes[k]){
//					printf("\n Boundary node %d is connected to another boundary node %d by %d bond ", boundary_nodes[i], boundary_nodes[k], j);
					occupation[boundary_nodes[i]*6+j]=1;
					occupation[connect[boundary_nodes[i]][j]*6+(j+3)%6]=1;					
				}
			}
		}
	}
 
	int bond_count_double = 0;   // these are bonds that will be double counted
	int bond_count_single = 0;   // these will be counted just once
	// counting how many bonds are there in the outer ring::
	for(i=0;i<outer_ring_count;i++){
		for (j=0;j<6;j++){
			dx_center = x[center] - x[connect[outer_ring_nodes[i]][j]];
			dy_center = y[center] - y[connect[outer_ring_nodes[i]][j]];

			dist_center = sqrt(dx_center*dx_center + dy_center*dy_center);

			if ((dist_center < (r_out-dr_center)) && (dist_center > r_in)){
				bond_count_double += 1;
			}
			else{
				bond_count_single += 1;
//				printf("\n single: %d ", connect[i][j]);
			}

			// finding the number of bonds in inner network for a p = 1 network
			if ((dist_center <= r_in)){
				bond_count_double += 1;
			}

		}
	}


	// counting how many bonds are there in the inner circle::
	int bond_count_inner_p1 = 0;
	int count_inner_bonds = 0;

	for(i=0;i<inner_count;i++){
		for (j=0;j<6;j++){
			dx_center = x[center] - x[connect[inner_nodes[i]][j]];
			dy_center = y[center] - y[connect[inner_nodes[i]][j]];
			dist_center = sqrt(dx_center*dx_center + dy_center*dy_center);

			// finding the number of bonds in inner network for a p = 1 network
			if ((dist_center <= r_in)){
				bond_count_inner_p1 += 1;

				if (occupation[i*6+j]==1){
					count_inner_bonds = count_inner_bonds+1;
				}

			}
		}
	}
	printf("\n For a p=1, network, there are %d total bonds in the inner circle \n", bond_count_inner_p1);
	printf("\n In this network, there are %d total bonds in the inner circle \n", count_inner_bonds);

	printf("\n For a p=1, network, there are %d double counted bonds in the outer ring \n", bond_count_double);
	printf("\n For a p=1, network, there are %d single counted bonds in the outer ring \n", bond_count_single);

	int outer_ring_bonds = (int)(bond_count_double/2. + bond_count_single);

	printf("\n For a p=1, network, there are %d total bonds in the outer ring \n", outer_ring_bonds);


//	for all simulations, pc flag was set to zero.
// counting total number of bonds
	int count_bonds = 0;

	// counting the bonds only in outer ring
	for(i=0;i<N;i++){

		dx_center = x[center] - x[i];
		dy_center = y[center] - y[i];
		dist_center = sqrt(dx_center*dx_center + dy_center*dy_center);

		if ((dist_center < (r_out-dr_center)) && (dist_center > r_in)){
			for(j=0;j<3;j++){
				if (occupation[i*6+j]==1){
					count_bonds = count_bonds+1;
				}
			}
		}
	}



//	int pc_num = (round) (PBOND*L*L*3);
	int pc_num = (round) (prob_eff*L*L*3);
	pc_num = (int)pc_num;
//	printf("\n Number of bonds required = %d PBOND = %f \n",pc_num,PBOND);
	printf("\n Number of bonds required = %d PBOND = %f \n",pc_num,PBOND);

	if(PBOND < 1.0){
		pc_num = (round) (PBOND*outer_ring_bonds);
		pc_num = (int)pc_num;

		printf("\n Total number of bonds int the outer ring is ::>> %d \n", count_bonds);
		printf("\n We need to remove %d bonds \n", (count_bonds-pc_num));
	}

	if ((count_bonds > pc_num) && (PBOND < 1.0))
	{
		printf("\n Total number of bonds is ::>> %d \n", count_bonds);
		printf("\n We need to remove %d bonds \n", (count_bonds-pc_num));
	}

//****************************************************************************************
// adding bonds to come to the right p value inside the outer ring only

// setting the coordination number to zero before counting
	for(i=0;i<N;i++){
		coord_num[i] = 0;
	}

// finding the coordination number again
	for(i=0;i<N;i++){
		for(j=0;j<6;j++){
			if (occupation[i*6+j]==1){
				coord_num[i] = coord_num[i]+1;
			}
		}
	}


	int rand_node, dip_node_flag, rand_bond, bond_flag;
	int temp_bond_count;

//****************************************************************************************

	// connecting the boundary nodes to each other::
	for(i=0;i<temp_count;i++){
		for (j=0;j<6;j++){
			for(k=0;k<temp_count;k++){
				if (connect[boundary_nodes[i]][j] == boundary_nodes[k]){
//					printf("\n Boundary node %d is connected to another boundary node %d by %d bond ", boundary_nodes[i], boundary_nodes[k], j);
					occupation[boundary_nodes[i]*6+j]=1;
					occupation[connect[boundary_nodes[i]][j]*6+(j+3)%6]=1;					
				}
			}
		}
	}

	count = 0;
	for(i=0;i<N;i++){
		for(j=0;j<3;j++){
			if (occupation[i*6+j]==1){
				count = count+1;
			}
		}
	}
	printf("\n Total number of bonds in the whole network is :::: %d \n", count);


// Final count of the coordination number
	for(i=0;i<N;i++){
		coord_num[i] = 0;
		for(j=0;j<6;j++){
			if (occupation[i*6+j]==1){
				coord_num[i] = coord_num[i]+1;
			}
		}
	}

	//***************************************************

//	int rand_node, dip_node_flag, kk, bond_flag, rand_bond;

	//realP=dangle();
	
	for(i=0;i<N;i++){
		for(j=0;j<6;j++){
			if( occupation[i*6+j]==1 && occupation[i*6+(j+1)%6]==1 ){
				if( rand_double() < PARP ){	
					arpoccupation[i*6+j]=1;
		
				}else{
					arpoccupation[i*6+j]=0;
				}
			}
		}
	}
//	dlattice(pr);
//***************************************************************
//************************SUB LATTICE****************************
//*********************FOR FIXED BOUNDRIES***********************

//Fixing the upper an lower edges removes 4*L degrees of freedom from the system.
		int rlen_flag = 0;

		times=10;  // times = times ; commented by abhinav for testing a smaller output file
//		times=1;


		rlen = 1.0;
		for(t=0;t<times;t++){
			if(t==0){

			output_boundary_nodes(boundary_nodes, temp_count);  // output circular boundary nodes
			output_inner_nodes(inner_nodes, inner_count);  // output circular boundary nodes			
			output_outer_nodes();  // output nodes outside the circular boundary and on the circular boundary
			output_dip_nodes();
			output_unique_dip_nodes();
			output_sub(pr);
			output_connect(pr);
			func_flag = 1;
			func(subpr);
			func_flag = 0;
			// calling a function to write all restlengths
			rlen_output();
			}

			rlen = rlen - ((1.0-rlen_max)/times);

			printf("rlen = %f \n", rlen);

			printf("\n steps %d mu %0.2f kappa2 %0.2e force %0.4f \n", times, mu, kappa2, force_val);
			printf("\n pbond %0.6f ranseed1 %ld ranseed2 %ld \n", PBOND, ranseed1, ranseed2);


//			setting up the rest length array
			for (i = 0; i < N; i++){
				for (j = 0; j < 6; j++){
					rlen_flag = 0;

					rlen_arr[i*6+j] = 1.0;
					// checking if we are at a dipole bond
					for (k=0;k<num_center;k++){
						row = (int)loc_center[k]/L;

						if (row%2 == 0){

							if ((i == loc_center[k]) || (connect[i][j] == loc_center[k])){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - 1) && (j == 2)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - 1) && (j == 4)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + 1) && (j == 1)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + 1) && (j == 5)){
								rlen_flag = 1;
							}

							else if ((i == loc_center[k] - L) && (j == 5)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - L) && (j == 3)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - L + 1) && (j == 0)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - L + 1) && (j == 4)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + L + 1) && (j == 2)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + L + 1) && (j == 0)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + L) && (j == 3)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + L) && (j == 1)){
								rlen_flag = 1;
							}
						
						}
						else{
							if ((i == loc_center[k]) || (connect[i][j] == loc_center[k])){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - 1) && (j == 2)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - 1) && (j == 4)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + 1) && (j == 1)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + 1) && (j == 5)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - L - 1) && (j == 5)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - L - 1) && (j == 3)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - L) && (j == 0)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] - L) && (j == 4)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + L) && (j == 2)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + L) && (j == 0)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + L - 1) && (j == 3)){
								rlen_flag = 1;
							}
							else if ((i == loc_center[k] + L - 1) && (j == 1)){
								rlen_flag = 1;
							}
						 }
					}

					if (rlen_flag == 1){
						rlen_arr[i*6+j] = rlen;
						rlen_arr[connect[i][j]*6 + (j+3)%6] = rlen;
					}
				}
			}

			// calling a function to write all restlengths
			rlen_output();

			// calling energy minimizer below ---------------------------

			frprmn(subpr,2*N,ftol,&iter,&fret,&func,&dfunc);
			func_flag = 1;
			ass=2*func(subpr)/N;
			func_flag = 0;

			// calling energy minimizer above ---------------------------

			// clamped:: ALL boundaries not being moved.
			for( i=0;i<N;i++ ){	

				dip_true = 0;    // flag for checking if we are at a dipole node

				for (jj=0;jj<(num_dip_unique);jj++){
					if ( (i == dip_node_unique[jj]) ) {
						dip_true = 1;
//						printf("\n we found a dipole!! :: %d \n",i);
					}
				}

		
//				if(i>=L && i<N-L && i%L != 0 && i%L != L-1 && dip_true == 0 ){								
//				if(i>=L && i<N-L && i%L != 0 && i%L != L-1){		// fixed rectangular boundary						
				if(outer_node[i] == 0){		// fixed circular boundary and everything outside it	
					pr[i]=subpr[i];
					pr[i+N]=subpr[i+N];
				}
			}

//**************************************************************
//*************WRITING OUTPUT***********************************			

//				dlattice(pr);
				output_sub(pr);
			j=0;
		}

		// writing the bending forces to file
		bend_force_output();
//***************************************************************
//*****************IMPLEMENTAION OF OTHER SUBROUTINES************
	
	}
}
//fclose(dq);

free(Ff);
free(dip_node);


}


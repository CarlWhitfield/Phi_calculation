//
//  FN_knot_code.h
//  
//
//  Created by Carl Whitfield on 17/05/2016.
//
//

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;

#define FROM_B_FILE 0
#define FROM_K_FILE 1
#define TORUS_KNOT 2
#define HOPF_LINK 3
#define BORR_RINGS 4
#define THREE_TWIST 5
#define FROM_UV_FILE 10
#define SIX_ONE 61
#define UNLINK 20

/*************************General maths and integer functions*****************************/

inline int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}

inline int incw(int i, int p, int N)    //increment with reflecting boundary between -1 and 0 and N-1 and N
{
    if(i+p<0) return (-(i+p+1));
    if(i+p>N-1) return (2*N-(i+p+1));
    return (i+p);
}

inline int sign(int i)
{
    if(i==0) return 0;
    else return i/abs(i);
}


/*************************Functions for knot initialisation*****************************/

double init_torus_knot(void);

double init_from_k_file(void);

double init_hopf_link(void);

double init_borr_rings(void);

double init_three_twist(void);

double init_six_one(void);

double init_unlink(void);

double initialise_knot();

//double add_twist(void);   //only works for a knot at the moment

/*************************Functions for B and Phi calcs*****************************/


 void initial_cond(double *x, double *y, double *z, double *phi, int *missed);
 
 void B_field_calc(double *x, double *y, double *z, double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed);
 
 void phi_calc(double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed, double *phi);
 
 int pathfind(int i0, int j0, int k0, int ie, int je, int ke, int *pi, int *pj, int *pk, int *ignore, double *Bx, double *By, double *Bz, double *Bmag);

//FitzHugh Nagumo functions
void uv_initialise(double *phi, double *u, double *v, double *ucv, int *missed);
void crossgrad_calc(double *u, double *v, double *ucv);
void uv_update(double *u, double *v);
void print_uv(double *u, double *v, double *ucv, double t);

/*************************File reading and writing*****************************/

int Bfile_read(double *phi, int *missed);
int uvfile_read(double *u,double *v);
void print_knot();
void print_B_phi(double *x, double *y, double*z, double *Bx, double *By, double *Bz, int *ignore, int *ignore1, int *missed, double *phi);
void print_uv(double *u, double *v, double *ucv, double t);
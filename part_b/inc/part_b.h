#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <time.h>

/* I/O directories */
// #define DATA_DIR '../data/'

/* physical constants in cgs units */
#define N_avo   6.0221367e23
#define k_boltz 1.380658e-16
#define a_rad   7.5646e-15
#define Msun    1.99e33
#define Rsun    6.96e10
#define Lsun    3.9e33
#define Tsun    5.78e3
#define G       6.67259e-8
#define c       2.99792458e10
#define sigma   5.67051e-5
/*
A Row of the opacity table. Contains opacity values for one value of log(T) and
several values of log(R).
*/
typedef struct{
    double * Rcol;
} OPAC_TROW;

void load_table(OPAC_TROW ** opac_Trow, double ** logT, double ** logR, int N_Trow,
    int N_Rcol);

double compute_mu_ion(double X, double Y, double Z);
double compute_mu_e(double X, double Y, double Z);
double compute_mu(double X, double Y, double Z);
double compute_density(double P, double T, double X, double Y, double Z);
double compute_opacity_escatt(double X);
double compute_opacity_ffabs(double X, double Y, double rho, double T);
double compute_opacity_bfabs(double X, double Z, double rho, double T);
double compute_opacity_H(double Z, double rho, double T);
double compute_opacity_approximate(double X, double Y, double Z, double rho, double T);
double compute_R(double T, double rho);
double compute_opacity(double X, double Y, double Z, double rho, double T,
    OPAC_TROW *opac_Trow, double *opac_logT, double *opac_logR, int N_Trow, int N_Rcol);
double pp_chain(double rho, double T, double X);
double cno_cycle(double rho, double T, double X, double Z);
double triple_alpha(double rho, double T, double Y);
double compute_energy(double rho, double T, double X, double Y, double Z);
double compute_P_gas(double P_tot, double T);
double compute_gamma2(double P, double T);
double compute_drdm(double rho, double r);
double compute_dPdm(double m, double r);
double compute_dLdm(double epsilon);
double compute_ratio(double kappa, double gamma2, double m, double L, double T, double P);
double compute_dTdm(double T, double P, double m, double r, double gamma2, double kappa, double L, double ratio);
double compute_center_radius(double m, double rho);
double compute_center_luminosity(double m, double epsilon);

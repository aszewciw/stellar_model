#include "part_a.h"

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* -------------------------------- MAIN --------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */


int main(int argc, char * argv[]){

    double logTc;   /* log10 of Central temperature */
    double logPc;   /* log10 of Central pressure */
    double Tc;      /* unlogged quantity */
    double Pc;
    int N_Rcol=19;  /* # of rows (logR values) in opacity table */
    int N_Trow=70;  /* # of columns (logT values) in opacity table */
    OPAC_TROW *opac_Trow;   /* table of log(opacity) values */
    double *opac_logT;      /* logT values of table */
    double *opac_logR;      /* logR values of table */
    double X;   /* Hydrogen mass fraction */
    double Y;   /* Helium mass fraction */
    double Z;   /* Metals mass fraction */
    double M;   /* Star mass */
    double P;   /* pressure */
    double T;   /* temperature */
    double rho; /* density */
    double kappa;   /* opacity */
    double epsilon; /* energy generation rate */
    double gamma2;  /* adiabatic term */
    double m;   /* enclosed mass */
    double r;   /* radius */
    double L;   /* luminosity */
    double dm;  /* mass step */
    double drdm; /* dr/dm */
    double dPdm; /* dP/dm */
    double dLdm; /* dL/dm */
    double dTdm; /* dT/dm */
    double dTdr_ratio;  /* {dT/dr}_rad / {dT/dr}_ad */
    int i;

    /* boundary conditions */
    logPc = 16.675;
    logTc = 7.46;
    Pc = pow(10.0, logPc);
    Tc = pow(10.0, logTc);

    /* Known star properties */
    X = 0.7;
    Y = 0.28;
    Z = 0.02;
    M = 7.08*Msun;

    /* choose step size */
    dm = M * pow(10.0,-5);
    /* initialize mass with nonzero value */
    m = dm;

    /* get opacity table */
    load_table(&opac_Trow, &opac_logT, &opac_logR, N_Trow, N_Rcol);

    /* compute initial values based on Pc and Tc */
    rho = compute_density(Pc, Tc, X, Y, Z);
    kappa = compute_opacity(X, Y, Z, rho, Tc, opac_Trow, opac_logT, opac_logR,
        N_Trow, N_Rcol);
    epsilon = compute_energy(rho, Tc, X, Y, Z);
    gamma2 = compute_gamma2(Pc, Tc);

    /* get r, L slightly off center and recompute values */
    r = compute_center_radius(m, rho);
    L = compute_center_luminosity(m, epsilon);

    /* get rad/ad ratio to see how to get temp */
    dTdr_ratio = compute_ratio(kappa, gamma2, m, L, Tc, Pc);
    P = Pc - (2.0/3.0)*M_PI*G*rho*rho*r*r;
    if(dTdr_ratio>1.0){ // convective
        T = Tc - (1.0-1.0/gamma2)*2.0*M_PI*G*rho*rho*Tc*r*r/Pc;
    }
    else{   // radiative
        T = Tc - (kappa*rho*rho*epsilon*r*r/(8.0*a_rad*c*pow(Tc,3.0)));
    }

    /* Recompute above quantities given adjustments */
    rho = compute_density(P, T, X, Y, Z);
    kappa = compute_opacity(X, Y, Z, rho, T, opac_Trow, opac_logT, opac_logR,
        N_Trow, N_Rcol);
    epsilon = compute_energy(rho, T, X, Y, Z);
    gamma2 = compute_gamma2(P, T);
    dTdr_ratio = compute_ratio(kappa, gamma2, m, L, T, P);

    /* File for output data */
    char filename[256];
    FILE *file;
    snprintf(filename, 256, "../data/part_a.dat");
    file = fopen(filename, "a");

    /* write column names and initial values */
    fprintf(file, "m(Msun)\tlogP\tlogT\tlogrho\tr(Rsun)\tL(Lsun)\topac\tenergy\tratio\n");
    fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
            m/Msun, log10(P), log10(T), log10(rho), r/Rsun, L/Lsun, kappa, epsilon, dTdr_ratio);

    /* integrate until we exceed the star's mass */
    while(m<M){

        /* compute all derivatives */
        drdm = compute_drdm(rho, r);
        dPdm = compute_dPdm(m, r);
        dLdm = compute_dLdm(epsilon);
        dTdr_ratio = compute_ratio(kappa, gamma2, m, L, T, P);
        dTdm = compute_dTdm(T, P, m, r, gamma2, kappa, L, dTdr_ratio);

        /* update values */
        r = r + drdm*dm;
        P = P + dPdm*dm;
        L = L + dLdm*dm;
        T = T + dTdm*dm;
        m = m + dm;

        /* calculate new quantities based on updated values */
        rho = compute_density(P,T,X,Y,Z);
        kappa = compute_opacity(X,Y,Z,rho,T,opac_Trow,opac_logT,opac_logR,N_Trow,N_Rcol);
        epsilon = compute_energy(rho,T,X,Y,Z);
        gamma2 = compute_gamma2(P,T);

        /* output step to file */
        fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                m/Msun, log10(P), log10(T), log10(rho), r/Rsun, L/Lsun, kappa, epsilon, dTdr_ratio);
    }
    fclose(file);

    /* Free allocated values */
    for(i=0; i<N_Trow; i++){
        free(opac_Trow[i].Rcol);
    }
    free(opac_Trow);
    free(opac_logR);
    free(opac_logT);

    return EXIT_SUCCESS;
}

/* ----------------------------------------------------------------------- */
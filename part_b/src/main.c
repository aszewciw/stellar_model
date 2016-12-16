#include "part_b.h"

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* -------------------------------- MAIN --------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */


int main(int argc, char * argv[]){

    /* first check that we pass the correct number of elements */
    if (argc != 6){
        fprintf( stderr, "Usage:\n %s log(Tc) log(Pc) log(Rs) log(Ls) write_flag > outfile\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    double logTc;   /* log10 of Central temperature */
    double logPc;   /* log10 of Central pressure */
    double logRs;   /* log10 of Star radius */
    double logLs;   /* log10 of Surface luminosity */
    double Tc;      /* unlogged quantity */
    double Pc;
    double Rs;
    double Ls;
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
    double mu;  /* mean molecular weight */
    int i;
    int write_flag;  /* 0 if we don't want to write data to file */

    /* get boundary conditions */
    logTc = atof(argv[1]);
    logPc = atof(argv[2]);
    logRs = atof(argv[3]);
    logLs = atof(argv[4]);
    write_flag = atoi(argv[5]);

    /* Convert to unlogged version */
    Tc = pow(10.0, logTc);
    Pc = pow(10.0, logPc);
    Rs = pow(10.0, logRs);
    Ls = pow(10.0, logLs);

    /* Known star properties */
    X = 0.7;
    Y = 0.28;
    Z = 0.02;
    M = 2.82*Msun;

    /* get mean molecular weight */
    mu = compute_mu(X,Y,Z);

    /* get opacity table */
    load_table(&opac_Trow, &opac_logT, &opac_logR, N_Trow, N_Rcol);

    /* compute initial quantities */
    rho = compute_density(Pc, Tc, X, Y, Z);
    kappa = compute_opacity(X, Y, Z, rho, Tc, opac_Trow, opac_logT, opac_logR,
        N_Trow, N_Rcol);
    epsilon = compute_energy(rho, Tc, X, Y, Z);
    gamma2 = compute_gamma2(Pc, Tc);

    /* choose step size */
    dm = M * pow(10.0,-5);
    /* initialize mass with nonzero value */
    m = dm;

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

    char filename[256];
    FILE *file;
    snprintf(filename, 256, "../data/forward.dat");
    file = fopen(filename, "a");

    if(write_flag==1){

        fprintf(file, "m(Msun)\tlogP\tlogT\tlogrho\tr(Rsun)\tL(Lsun)\topac\tenergy\tratio\n");
        fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                m/Msun, log10(P), log10(T), log10(rho), r/Rsun, L/Lsun, kappa, epsilon, dTdr_ratio);
    }

    /* Perform forward integration */
    while(m<M/2.0){
        drdm = compute_drdm(rho, r);
        dPdm = compute_dPdm(m, r);
        dLdm = compute_dLdm(epsilon);
        dTdr_ratio = compute_ratio(kappa, gamma2, m, L, T, P);
        dTdm = compute_dTdm(T, P, m, r, gamma2, kappa, L, dTdr_ratio);

        r = r + drdm*dm;
        P = P + dPdm*dm;
        L = L + dLdm*dm;
        T = T + dTdm*dm;
        m = m + dm;

        rho = compute_density(P,T,X,Y,Z);
        kappa = compute_opacity(X,Y,Z,rho,T,opac_Trow,opac_logT,opac_logR,N_Trow,N_Rcol);
        epsilon = compute_energy(rho,T,X,Y,Z);
        gamma2 = compute_gamma2(P,T);
        if(write_flag==1){
            fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                    m/Msun, log10(P), log10(T), log10(rho), r/Rsun, L/Lsun, kappa, epsilon, dTdr_ratio);
        }

    }
    fclose(file);


    double logP_low, logT_low, logr_low, logL_low;
    logP_low = log10(P);
    logT_low = log10(T);
    logr_low = log10(r);
    logL_low = log10(L);

    /* now reset things to the surface and integrate backwards */
    m = M;
    r = Rs;
    L = Ls;
    T = pow(L/(4.0*M_PI*r*r*sigma), 0.25);
    rho = 1.0e-5;
    P = a_rad*pow(T,4.0)/3.0 + N_avo*k_boltz*rho*T/mu;
    // rho = compute_density(P,T,X,Y,Z);
    kappa = compute_opacity(X,Y,Z,rho,T,opac_Trow,opac_logT,opac_logR,N_Trow,N_Rcol);
    epsilon = compute_energy(rho,T,X,Y,Z);
    gamma2 = compute_gamma2(P,T);
    dTdr_ratio = compute_ratio(kappa, gamma2, m, L, T, P);

    snprintf(filename, 256, "../data/backward.dat");
    file = fopen(filename, "a");

    if(write_flag==1){
        fprintf(file, "m(Msun)\tlogP\tlogT\tlogrho\tr(Rsun)\tL(Lsun)\topac\tenergy\tratio\n");
        fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                m/Msun, log10(P), log10(T), log10(rho), r/Rsun, L/Lsun, kappa, epsilon, dTdr_ratio);
    }

    /* Perform backward integration */
    while(m>M/2.0){
        drdm = compute_drdm(rho, r);
        dPdm = compute_dPdm(m, r);
        dLdm = compute_dLdm(epsilon);
        dTdr_ratio = compute_ratio(kappa, gamma2, m, L, T, P);
        dTdm = compute_dTdm(T, P, m, r, gamma2, kappa, L, dTdr_ratio);

        r = r - drdm*dm;
        P = P - dPdm*dm;
        L = L - dLdm*dm;
        T = T - dTdm*dm;
        m = m - dm;

        rho = compute_density(P,T,X,Y,Z);
        kappa = compute_opacity(X,Y,Z,rho,T,opac_Trow,opac_logT,opac_logR,N_Trow,N_Rcol);
        epsilon = compute_energy(rho,T,X,Y,Z);
        gamma2 = compute_gamma2(P,T);
        if(write_flag==1){
            fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                    m/Msun, log10(P), log10(T), log10(rho), r/Rsun, L/Lsun, kappa, epsilon, dTdr_ratio);
        }

    }
    fclose(file);

    double logP_high, logT_high, logr_high, logL_high;
    logP_high = log10(P);
    logT_high = log10(T);
    logr_high = log10(r);
    logL_high = log10(L);

    double dlogP = logP_high - logP_low;
    double dlogT = logT_high - logT_low;
    double dlogr = logr_high - logr_low;
    double dlogL = logL_high - logL_low;

    fprintf(stdout, "%lf\t%lf\t%lf\t%lf\n", dlogP, dlogT, dlogr, dlogL);

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
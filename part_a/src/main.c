#include "part_a.h"

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* -------------------------------- MAIN --------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */


int main(int argc, char * argv[]){

    int N_Rcol = 19;
    int N_Trow = 70;
    OPAC_TROW *opac_Trow;
    double *opac_logT;
    double *opac_logR;
    double X,Y,Z,M,P,T;
    int i;
    double logTc, logPc;

    load_table(&opac_Trow, &opac_logT, &opac_logR, N_Trow, N_Rcol);

    X = 0.7;
    Y = 0.28;
    Z = 0.02;

    logPc = 16.675;
    logTc = 7.46;
    M = 7.08*Msun;
    double Pc = pow(10.0, logPc);
    double Tc = pow(10.0, logTc);

    double rho = compute_density(Pc, Tc, X, Y, Z);
    fprintf(stderr, "Starting density is %lf g/cm^3\n", rho);

    double kappa = compute_opacity(X, Y, Z, rho, Tc, opac_Trow, opac_logT, opac_logR,
        N_Trow, N_Rcol);
    fprintf(stderr, "Starting opacity is %lf cm^2/g\n", kappa);

    double epsilon = compute_energy(rho, Tc, X, Y, Z);
    fprintf(stderr, "Starting energy is %lf ergs\n", epsilon);

    double gamma2 = compute_gamma2(Pc, Tc);
    fprintf(stderr, "Starting gamma2 is %lf\n", gamma2);

    double m, r, L, dm;
    double drdm, dPdm, dLdm, dTdm, dTdr_ratio;
    dm = M * pow(10.0,-5);

    // m = dm / 1.0e5;
    m = dm;

    r = compute_center_radius(m, rho);
    L = compute_center_luminosity(m, epsilon);

    // I think maybe we need to do this adjustment?
    P = Pc - (2.0/3.0)*M_PI*G*rho*rho*r*r;
    T = Tc - (1.0-1.0/gamma2)*2.0*M_PI*G*rho*rho*Tc*r*r/Pc;

    rho = compute_density(P, T, X, Y, Z);
    kappa = compute_opacity(X, Y, Z, rho, T, opac_Trow, opac_logT, opac_logR,
        N_Trow, N_Rcol);
    epsilon = compute_energy(rho, T, X, Y, Z);
    gamma2 = compute_gamma2(P, T);
    dTdr_ratio = compute_ratio(kappa, gamma2, m, L, T, P);


    char filename[256];
    FILE *file;
    snprintf(filename, 256, "../data/output.dat");
    file = fopen(filename, "a");
    fprintf(file, "m(Msun)\tlogP\tlogT\tlogrho\tr(Rsun)\tL(Lsun)\topac\tenergy\tratio\n");
    fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
            m/Msun, log10(P), log10(T), log10(rho), r/Rsun, L/Lsun, kappa, epsilon, dTdr_ratio);

    while(m<M){
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
        fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                m/Msun, log10(P), log10(T), log10(rho), r/Rsun, L/Lsun, kappa, epsilon, dTdr_ratio);
    }
    fclose(file);

    fprintf(stderr, "Printing final quantities in solar units.\n");
    fprintf(stderr, "Final mass: %lf\n", m/Msun);
    fprintf(stderr, "Final luminosity: %lf\n", L/Lsun);
    fprintf(stderr, "Final radius: %lf\n", r/Rsun);
    fprintf(stderr, "Final log(T): %lf\n", log10(T));


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
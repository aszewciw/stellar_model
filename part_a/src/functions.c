#include "part_a.h"

/* ----------------------------------------------------------------------- */
/* --------------------  Density-related functions  ---------------------- */
/* ----------------------------------------------------------------------- */

/* compute mu_ion */
double compute_mu_ion(double X, double Y, double Z){
    /*
    X, Y, Z - mass fractions of H, He, metals

    Currently we make the assumption that Z / A_Z = 0
    */

    double mu;         /* ion mean molecular weight */
    double A_X, A_Y;   /* nuclear mass numbers */

    A_X = 1.0;
    A_Y = 4.0;

    mu = 1.0 / ( (X / A_X) + (Y / A_Y) );

    return mu;
}
/* ----------------------------------------------------------------------- */

/* compute mu_e */
double compute_mu_e(double X, double Y, double Z){
    /*
    X, Y, Z -       mass fractions of H, He, metals
    y_X, y_Y, y_Z - ionization fractions

    Assume that A_Z / Q_Z = 2.
    */

    double mu;          /* ion mean molecular weight */
    double A_X, A_Y;    /* nuclear mass numbers */
    double Q_X, Q_Y;    /* nuclear charge numbers */

    A_X = 1.0;
    A_Y = 4.0;
    Q_X = 1.0;
    Q_Y = 2.0;

    mu = 1.0 / ( (X * Q_X / A_X) + (Y * Q_Y / A_Y) + (Z * 0.5) );

    return mu;
}
/* ----------------------------------------------------------------------- */

/* compute mean molecular weight */
double compute_mu(double X, double Y, double Z){
    /*
    X, Y, Z -       mass fractions of H, He, metals
    y_X, y_Y, y_Z - ionization fractions
    */

    double mu, mu_i, mu_e;

    mu_i = compute_mu_ion(X, Y, Z);
    mu_e = compute_mu_e(X, Y, Z);

    mu = 1.0 / ( (1.0/mu_i) + (1.0/mu_e) );

    return mu;
}
/* ----------------------------------------------------------------------- */

double compute_density(double P, double T, double X, double Y, double Z){

    /* make sure pressure doesn't drop below (1/3)*a*T^4 */

    double mu = compute_mu(X, Y, Z);
    double rho;

    rho = (P - a_rad*pow(T,4.0)/3.0) * (mu/(N_avo*k_boltz*T));

    return rho;
}
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/* --------------------  Opacity-related functions  ---------------------- */
/* ----------------------------------------------------------------------- */

double BilinearInterpolation(double q11, double q12, double q21, double q22,
    double x1, double x2, double y1, double y2, double x, double y)
{

    double x2x1, y2y1, x2x, y2y, yy1, xx1;

    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;

    return 1.0 / (x2x1 * y2y1) * (
        q11 * x2x * y2y +
        q21 * xx1 * y2y +
        q12 * x2x * yy1 +
        q22 * xx1 * yy1
    );
}

/* ----------------------------------------------------------------------- */

double compute_opacity_escatt(double X){
    return 0.2*(1.0+X);
}
/* ----------------------------------------------------------------------- */

double compute_opacity_ffabs(double X, double Y, double rho, double T){
    return 4.0e22 * (X+Y) * (1.0+X) * rho * pow(T, -3.5);
}
/* ----------------------------------------------------------------------- */

double compute_opacity_bfabs(double X, double Z, double rho, double T){
    return 4.0e25 * Z * (1.0+X) * rho * pow(T, -3.5);
}
/* ----------------------------------------------------------------------- */

double compute_opacity_H(double Z, double rho, double T){
    return 2.5e-31 * (Z/0.02) * pow(rho, 0.5) * pow(T,9);
}
/* ----------------------------------------------------------------------- */

double compute_opacity_approximate(double X, double Y, double Z, double rho,
    double T){

    double k_e, k_ff, k_bf, k_H, k_tot, k_tmp;

    k_e = compute_opacity_escatt(X);
    k_ff = compute_opacity_ffabs(X,Y,rho,T);
    k_bf = compute_opacity_bfabs(X,Z,rho,T);
    k_H = compute_opacity_H(Z,rho,T);

    k_tmp = k_e + k_ff + k_bf;
    k_tot = 1.0/(1.0/k_H + 1.0/k_tmp);

    return k_tot;
}

/* ----------------------------------------------------------------------- */

double compute_R(double T, double rho){

    double R;

    R = rho / pow((T/1.0e6),3);

    return R;
}

/* ----------------------------------------------------------------------- */

double compute_opacity(double X, double Y, double Z, double rho, double T,
    OPAC_TROW *opac_Trow, double *opac_logT, double *opac_logR, int N_Trow, int N_Rcol){

    double R, logT, logR, kappa;
    double T1, T2, R1, R2, k11, k12, k21, k22;
    int i, j;
    int flag = 0;

    R = compute_R(T, rho);
    logT = log10(T);
    logR = log10(R);

    i = 0;
    while(i<N_Trow-1){
        T1 = opac_logT[i];
        T2 = opac_logT[i+1];
        if((logT>T1)&&(logT<T2)){
            break;
        }
        i++;
    }

    j = 0;
    while(j<N_Rcol-1){
        R1 = opac_logR[j];
        R2 = opac_logR[j+1];
        if( (logR>R1)&&(logR<R2) ){
            break;
        }
        j++;
    }

    if( (i==N_Trow-1)||(j==N_Rcol-1) ){
        flag = 1;
    }
    else{
        k11 = opac_Trow[i].Rcol[j];
        k12 = opac_Trow[i].Rcol[j+1];
        k21 = opac_Trow[i+1].Rcol[j];
        k22 = opac_Trow[i+1].Rcol[j+1];
        if( (k11==-99999.0) || (k12==-99999.0) || (k21==-99999.0) || (k22==-99999.0) ){
            flag = 1;
        }
    }

    if(flag==1){
        kappa = compute_opacity_approximate(X, Y, Z, rho, T);
    }
    else{
        kappa = BilinearInterpolation(k11, k12, k21, k22, T1, T2, R1, R2, logT, logR);
        kappa = pow(10.0,kappa);
    }

    return kappa;

}
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/* ---------------------  Energy-related functions  ---------------------- */
/* ----------------------------------------------------------------------- */

double pp_chain(double rho, double T, double X){

    double T9 = T/1.0e9;
    double energy;

    energy = ( 2.4e4 * rho * X * X * exp( -3.38 / pow(T9,(1.0/3.0)) )
                / pow(T9,(2.0/3.0)) );
    return energy;
}

/* ----------------------------------------------------------------------- */

double cno_cycle(double rho, double T, double X, double Z){

    double T9 = T/1.0e9;
    double energy;

    energy = ( 4.4e25 * rho * X * Z * exp( -15.228 / pow(T9,(1.0/3.0)) )
                / pow(T9,(2.0/3.0)) );
    return energy;
}
/* ----------------------------------------------------------------------- */

double triple_alpha(double rho, double T, double Y){

    double T9 = T/1.0e9;
    double energy;

    energy = ( 5.0e8 * rho * rho * Y * Y * Y * exp(-4.4 / T9)
                / pow(T9,(3.0)) );

    return energy;
}
/* ----------------------------------------------------------------------- */

double compute_energy(double rho, double T, double X, double Y, double Z){

    double e_tot, e_pp, e_cno, e_alpha;

    e_pp = pp_chain(rho, T, X);
    e_cno = cno_cycle(rho, T, X, Z);
    e_alpha = triple_alpha(rho, T, Y);

    e_tot = e_pp + e_cno + e_alpha;

    return e_tot;
}
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* ----------------------  Gamma-related functions  ---------------------- */
/* ----------------------------------------------------------------------- */

double compute_P_gas(double P_tot, double T){

    /* Get gas pressure by subtracting radiation pressure */

    double P_rad, P_gas;

    P_rad = a_rad * pow(T, 4.0) / 3.0;

    P_gas = P_tot - P_rad;

    if(P_gas < 0.0){
        fprintf(stderr, "Error! Gas pressure is negative! Exiting code...\n");
        exit(-1);
    }

    return P_gas;
}
/* ----------------------------------------------------------------------- */

double compute_gamma2(double P, double T){

    double B;
    double P_gas;
    double gamma2;

    P_gas = compute_P_gas(P, T);
    B = P_gas / P;

    gamma2 = (32.0 - 24.0*B - 3.0*B*B) / (24.0 - 18.0*B - 3.0*B*B);

    return gamma2;
}

/* ----------------------------------------------------------------------- */

double compute_drdm(double rho, double r){
    return(0.25 / (M_PI * r * r * rho));
}

double compute_dPdm(double m, double r){
    return( (-0.25*G*m)/(M_PI*pow(r,4)) );
}

double compute_dLdm(double epsilon){
    return epsilon;
}


double compute_ratio(double kappa, double gamma2, double m, double L, double T, double P){

    double ratio;

    ratio = 3.0*kappa*L*P/(1.0-1.0/gamma2)/16.0/M_PI/a_rad/c/pow(T,4)/G/m;

    return ratio;
}

double compute_dTdm(double T, double P, double m, double r, double gamma2, double kappa, double L, double ratio){

    double dTdm;

    /* check convection condition */
    if(ratio>1.0){
        dTdm = -(1.0 - 1.0/gamma2)*G*m*T/(4.0*M_PI*pow(r,4)*P);
    }
    else{
        dTdm = -3.0*kappa*L/(64.0*M_PI*M_PI*a_rad*c*pow(r,4)*pow(T,3));
    }

    return dTdm;
}

double compute_center_radius(double m, double rho){
    return pow( (0.75*m/M_PI/rho), (1.0/3.0) );
}

double compute_center_luminosity(double m, double epsilon){
    return epsilon*m;
}
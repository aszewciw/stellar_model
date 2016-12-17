#include "part_a.h"

/* load opacity table */
void load_table(OPAC_TROW ** opac_Trow, double ** logT, double ** logR, int N_Trow,
    int N_Rcol)
{
    int i, j;
    char opac_filename[256];
    FILE *opac_file;
    OPAC_TROW *o;
    double *R;
    double *T;
    double *kappa;

    snprintf(opac_filename, 256, "../data/OPAL_X0.7_Y0.28.dat");

    if((opac_file=fopen(opac_filename,"r"))==NULL){
        fprintf(stderr, "Error: Cannot open file %s \n", opac_filename);
        exit(EXIT_FAILURE);
    }

    /* Claim array for  */
    o = calloc(N_Trow, sizeof(OPAC_TROW));
    R = calloc(N_Rcol, sizeof(double));
    T = calloc(N_Trow, sizeof(double));

    i = 0;
    /* Read in first row of log(R) values */
    while(i<N_Rcol){
        double R_tmp;
        fscanf(opac_file, "%lf", &R_tmp);

        /* skip the first element */
        if(R_tmp==9999.0) continue;

        R[i] = R_tmp;
        i++;
    }

    for(i=0; i<N_Trow; i++){

        fscanf(opac_file, "%lf", &T[i]);
        kappa = calloc(N_Rcol, sizeof(double));

        for(j=0; j<N_Rcol; j++){
            fscanf(opac_file, "%lf", &kappa[j]);
        }

        o[i].Rcol = kappa;
    }
    fclose(opac_file);

    /* Assign values to main function arguments */
    *logT = T;
    *logR = R;
    *opac_Trow = o;
}


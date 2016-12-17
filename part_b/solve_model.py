# import pandas as pd
import os, sys
# import matplotlib.pyplot as plt
import numpy as np

# Get partial derivative
def partialY_partialX(Y_plus, Y_minus, dX):

    return (Y_plus - Y_minus)/(2.0*dX)

def main():

    # starting values
    logTc =  7.36005299771
    logPc =  17.0933196124
    logRs =  11.2068629231
    logLs =  35.3504487033

    # control how we calculate the derivatives
    fraction = 0.001

    # step sizes
    dlogT = fraction*logTc
    dlogP = fraction*logPc
    dlogR = fraction*logRs
    dlogL = fraction*logLs

    # control size of step for new test values (after computing derivatives)
    step_fraction = 0.1

    # tolerances for ending loop
    log_tol = 0.001

    '''
    See the main.c code for more details.
    What is saved in these files are the differences between forward integration
    and backward integration at the center for the following quantites:
        log(T), log(P), log(R), log(L)
    No reference is made in the file to this, so the following code is not
    elegant. It exploits that this is the data structure.
    '''

    i = 0   # step counter
    step_limit = 100    # set limit to avoid infinite loop

    while(i<step_limit):

        # print info
        print('On step ' + str(i))
        print("logTc = ", logTc)
        print("logPc = ", logPc)
        print("logRs = ", logRs)
        print("logLs = ", logLs)

        # perform trial; called "fiducial"
        fname = "../data/trial.dat"
        cmd = ( "./bin/part_b " + str(logTc) + ' ' + str(logPc) + ' ' + str(logRs)
                + ' ' + str(logLs) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        fiducial = np.genfromtxt(fname)

        print("delta(LogTc) = ", fiducial[0])
        print("delta(LogPc) = ", fiducial[1])
        print("delta(LogRs) = ", fiducial[2])
        print("delta(LogLs) = ", fiducial[3])

        # Check if we're within our tolerance threshhold and end if so
        if(np.all(np.fabs(fiducial)<log_tol)):
            fname = "../data/trial.dat"
            cmd = ( "./bin/part_b " + str(logTc) + ' ' + str(logPc) + ' ' + str(logRs)
                    + ' ' + str(logLs) + ' ' + str(1) + ' > ' + fname )
            os.system(cmd)
            break

        # Calculate new trial values
        logT_m = logTc - dlogT
        logT_p = logTc + dlogT
        logP_m = logPc - dlogP
        logP_p = logPc + dlogP
        logR_m = logRs - dlogR
        logR_p = logRs + dlogR
        logL_m = logLs - dlogL
        logL_p = logLs + dlogL

        # New steps in T
        fname = "../data/Tminus.dat"
        cmd = ( "./bin/part_b " + str(logT_m) + ' ' + str(logPc) + ' ' + str(logRs)
                + ' ' + str(logLs) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        Tminus = np.genfromtxt(fname)

        fname = "../data/Tplus.dat"
        cmd = ( "./bin/part_b " + str(logT_p) + ' ' + str(logPc) + ' ' + str(logRs)
                + ' ' + str(logLs) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        Tplus = np.genfromtxt(fname)

        # New steps in P
        fname = "../data/Pminus.dat"
        cmd = ( "./bin/part_b " + str(logTc) + ' ' + str(logP_m) + ' ' + str(logRs)
                + ' ' + str(logLs) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        Pminus = np.genfromtxt(fname)

        fname = "../data/Pplus.dat"
        cmd = ( "./bin/part_b " + str(logTc) + ' ' + str(logP_p) + ' ' + str(logRs)
                + ' ' + str(logLs) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        Pplus = np.genfromtxt(fname)

        # New steps in R
        fname = "../data/Rminus.dat"
        cmd = ( "./bin/part_b " + str(logTc) + ' ' + str(logPc) + ' ' + str(logR_m)
                + ' ' + str(logLs) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        Rminus = np.genfromtxt(fname)

        fname = "../data/Rplus.dat"
        cmd = ( "./bin/part_b " + str(logTc) + ' ' + str(logPc) + ' ' + str(logR_p)
                + ' ' + str(logLs) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        Rplus = np.genfromtxt(fname)

        # New steps in L
        fname = "../data/Lminus.dat"
        cmd = ( "./bin/part_b " + str(logTc) + ' ' + str(logPc) + ' ' + str(logRs)
                + ' ' + str(logL_m) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        Lminus = np.genfromtxt(fname)

        fname = "../data/Lplus.dat"
        cmd = ( "./bin/part_b " + str(logTc) + ' ' + str(logPc) + ' ' + str(logRs)
                + ' ' + str(logL_p) + ' ' + str(0) + ' > ' + fname )
        os.system(cmd)
        Lplus = np.genfromtxt(fname)

        # Difference matrix
        # Ordering of rows and columns is T,P,R,L
        diff_mat = np.zeros((4,4))
        diff_mat[0,0] = partialY_partialX(Tplus[0], Tminus[0], dlogT)
        diff_mat[0,1] = partialY_partialX(Tplus[1], Tminus[1], dlogP)
        diff_mat[0,2] = partialY_partialX(Tplus[2], Tminus[2], dlogR)
        diff_mat[0,3] = partialY_partialX(Tplus[3], Tminus[3], dlogL)

        diff_mat[1,0] = partialY_partialX(Pplus[0], Pminus[0], dlogT)
        diff_mat[1,1] = partialY_partialX(Pplus[1], Pminus[1], dlogP)
        diff_mat[1,2] = partialY_partialX(Pplus[2], Pminus[2], dlogR)
        diff_mat[1,3] = partialY_partialX(Pplus[3], Pminus[3], dlogL)

        diff_mat[2,0] = partialY_partialX(Rplus[0], Rminus[0], dlogT)
        diff_mat[2,1] = partialY_partialX(Rplus[1], Rminus[1], dlogP)
        diff_mat[2,2] = partialY_partialX(Rplus[2], Rminus[2], dlogR)
        diff_mat[2,3] = partialY_partialX(Rplus[3], Rminus[3], dlogL)

        diff_mat[3,0] = partialY_partialX(Lplus[0], Lminus[0], dlogT)
        diff_mat[3,1] = partialY_partialX(Lplus[1], Lminus[1], dlogP)
        diff_mat[3,2] = partialY_partialX(Lplus[2], Lminus[2], dlogR)
        diff_mat[3,3] = partialY_partialX(Lplus[3], Lminus[3], dlogL)

        # Determine new step by inverting matrix
        inv_mat = np.linalg.inv(diff_mat)
        steps = inv_mat.dot(-fiducial)
        steps *= step_fraction

        # Set new values
        logTc+=steps[0]
        logPc+=steps[1]
        logRs+=steps[2]
        logLs+=steps[3]

        i = i+1

    # Plot now that we've converged
    cmd = 'python plot_results.py'
    os.system(cmd)

if __name__ == '__main__':
    main()
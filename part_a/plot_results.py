import pandas as pd
import matplotlib.pyplot as plt

def main():

    filename = '../data/output.dat'
    plots_dir = './plots/'

    df = pd.read_csv(filename, sep='\s+')

    plt.clf()
    plt.figure(1)
    plt.plot(df['m(Msun)'], df['logP'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('log(P)')
    plt.savefig(plots_dir + 'P_vs_m.png')

    plt.clf()
    plt.figure(2)
    plt.plot(df['m(Msun)'], df['logT'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('log(T)')
    plt.savefig(plots_dir + 'T_vs_m.png')

    plt.clf()
    plt.figure(3)
    plt.plot(df['m(Msun)'], df['logrho'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('log(rho)')
    plt.savefig(plots_dir + 'rho_vs_m.png')

    plt.clf()
    plt.figure(4)
    plt.plot(df['m(Msun)'], df['r(Rsun)'])
    # plt.plot(df['r(Rsun)'], df['m(Msun)'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('Radius (Rsun)')
    plt.savefig(plots_dir + 'r_vs_m.png')

    plt.clf()
    plt.figure(5)
    plt.plot(df['m(Msun)'], df['L(Lsun)'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('L (Lsun)')
    plt.savefig(plots_dir + 'L_vs_m.png')

    plt.clf()
    plt.figure(6)
    plt.plot(df['m(Msun)'], df['opac'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('kappa')
    plt.axis([0, 7, 0, 10])
    plt.savefig(plots_dir + 'opac_vs_m.png')

    plt.clf()
    plt.figure(7)
    plt.plot(df['m(Msun)'], df['energy'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('energy rate')
    plt.savefig(plots_dir + 'E_vs_m.png')

    plt.clf()
    plt.figure(1)
    plt.plot(df['m(Msun)'], df['ratio'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('ratio')
    plt.axis([0, 7, 0, 10])
    plt.savefig(plots_dir + 'ratio_vs_m.png')

if __name__ == '__main__':
    main()
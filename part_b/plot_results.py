import pandas as pd
import matplotlib.pyplot as plt

def main():

    plots_dir = './plots/'

    filename = '../data/forward.dat'
    df_for = pd.read_csv(filename, sep='\s+')

    filename = '../data/backward.dat'
    df_back = pd.read_csv(filename, sep='\s+')

    plt.clf()
    plt.figure(1)
    plt.plot(df_for['m(Msun)'], df_for['logP'])
    plt.plot(df_back['m(Msun)'], df_back['logP'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('log(P)')
    plt.savefig(plots_dir + 'P_vs_m.png')

    plt.clf()
    plt.figure(2)
    plt.plot(df_for['m(Msun)'], df_for['logT'])
    plt.plot(df_back['m(Msun)'], df_back['logT'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('log(T)')
    plt.savefig(plots_dir + 'T_vs_m.png')

    plt.clf()
    plt.figure(3)
    plt.plot(df_for['m(Msun)'], df_for['logrho'])
    plt.plot(df_back['m(Msun)'], df_back['logrho'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('log(rho)')
    plt.savefig(plots_dir + 'rho_vs_m.png')

    plt.clf()
    plt.figure(4)
    plt.plot(df_for['m(Msun)'], df_for['r(Rsun)'])
    plt.plot(df_back['m(Msun)'], df_back['r(Rsun)'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('Radius (Rsun)')
    plt.savefig(plots_dir + 'r_vs_m.png')

    plt.clf()
    plt.figure(5)
    plt.plot(df_for['m(Msun)'], df_for['L(Lsun)'])
    plt.plot(df_back['m(Msun)'], df_back['L(Lsun)'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('L (Lsun)')
    plt.savefig(plots_dir + 'L_vs_m.png')

    plt.clf()
    plt.figure(6)
    plt.plot(df_for['m(Msun)'], df_for['opac'])
    plt.plot(df_back['m(Msun)'], df_back['opac'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('kappa')
    plt.axis([0, 3, 0, 10])
    plt.savefig(plots_dir + 'opac_vs_m.png')

    plt.clf()
    plt.figure(7)
    plt.plot(df_for['m(Msun)'], df_for['energy'])
    plt.plot(df_back['m(Msun)'], df_back['energy'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('energy rate')
    plt.savefig(plots_dir + 'E_vs_m.png')

    plt.clf()
    plt.figure(8)
    plt.plot(df_for['m(Msun)'], df_for['ratio'])
    plt.plot(df_back['m(Msun)'], df_back['ratio'])
    plt.xlabel('Mass (Msun)')
    plt.ylabel('ratio')
    plt.axis([0, 3, 0, 10])
    plt.savefig(plots_dir + 'ratio_vs_m.png')

if __name__ == '__main__':
    main()
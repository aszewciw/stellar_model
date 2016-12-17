import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    '''
    Produce all plots from part a
    '''

    filename = '../data/part_a.dat'
    plots_dir = './plots/'

    df = pd.read_csv(filename, sep='\s+')

    plt.clf()
    plt.figure(1)
    plt.plot(df['m(Msun)'], df['logP'])
    plt.xlabel(r'$\frac{m}{M_{sun}}$', fontsize=20)
    plt.ylabel(r'$log(P)$, (Ba)')
    plt.title('Pressure vs. mass')
    plt.savefig(plots_dir + 'P_vs_m.png')

    plt.clf()
    plt.figure(2)
    plt.plot(df['m(Msun)'], df['logT'])
    plt.xlabel(r'$\frac{m}{M_{sun}}$', fontsize=20)
    plt.ylabel(r'$log(T)$, (K)')
    plt.title('Temperature vs. mass')
    plt.savefig(plots_dir + 'T_vs_m.png')

    plt.clf()
    plt.figure(3)
    plt.plot(df['m(Msun)'], df['logrho'])
    plt.xlabel(r'$\frac{m}{M_{sun}}$', fontsize=20)
    plt.ylabel(r'$log(\rho), (g/cm^3)$')
    plt.title('Density vs. mass')
    plt.savefig(plots_dir + 'rho_vs_m.png')

    plt.clf()
    plt.figure(4)
    plt.plot(df['m(Msun)'], df['r(Rsun)'])
    plt.xlabel(r'$\frac{m}{M_{sun}}$', fontsize=20)
    plt.ylabel(r'$\frac{r}{R_{sun}}$', fontsize=20)
    plt.title('Radius vs. mass')
    plt.savefig(plots_dir + 'r_vs_m.png')

    plt.clf()
    plt.figure(5)
    plt.plot(df['m(Msun)'], df['L(Lsun)'])
    plt.xlabel(r'$\frac{m}{M_{sun}}$', fontsize=20)
    plt.ylabel(r'$\frac{L}{L_{sun}}$', fontsize=20)
    plt.title('Luminosity vs. mass')
    plt.savefig(plots_dir + 'L_vs_m.png')

    plt.clf()
    plt.figure(6)
    plt.plot(df['m(Msun)'], df['opac'])
    plt.xlabel(r'$\frac{m}{M_{sun}}$', fontsize=20)
    plt.ylabel(r'$\kappa, cm^2/g$')
    plt.title('Opacity vs. mass')
    plt.axis([0, 7, 0, 10])
    plt.savefig(plots_dir + 'opac_vs_m.png')

    plt.clf()
    plt.figure(7)
    plt.plot(df['m(Msun)'], df['energy'])
    plt.xlabel(r'$\frac{m}{M_{sun}}$', fontsize=20)
    plt.ylabel('energy generation rate (ergs/s/g)')
    plt.title('Energy generation rate vs. mass')
    plt.savefig(plots_dir + 'E_vs_m.png')

    x = np.linspace(0, 7, 1000)
    y = np.ones(1000)
    plt.clf()
    plt.figure(8)
    plt.plot(df['m(Msun)'], df['ratio'])
    plt.plot(x,y,color='red')
    plt.xlabel(r'$\frac{m}{M_{sun}}$', fontsize=20)
    plt.ylabel(r'$(\frac{dT}{dr})_{rad}/(\frac{dT}{dr})_{ad}$', fontsize=24)
    plt.title('dT/dr ratio vs. mass')
    plt.axis([0, 7, 0, 10])
    plt.savefig(plots_dir + 'ratio_vs_m.png')

if __name__ == '__main__':
    main()
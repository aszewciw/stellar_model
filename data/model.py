# Stellar Model part 1
import numpy as np
from matplotlib import pyplot as plt
import math
from scipy import interpolate

# Constants
G = 6.67*10**(-8)
a = 7.5*10**(-15)
c = 3.*10**10
NA = 6.02*10**23
h = 6.63*10**(-27)
me = 9.1*10**(-28)
mp = 1.67*10**(-24)
k = 1.38*10**(-16)
sigma = 5.67*10**(-5)
Msun = 1.99*10**33
Lsun = 3.9*10**33
Rsun = 6.96*10**10

# Conditions
M = 7.08*Msun
X = 0.7
Y = 0.28
Z = 1. - (X + Y)
dm = M*10**(-5.)	# Step size
mu = (((X/1.) + (Y/4.)) + (((1.*X*1.)/1.)+((2.*Y*1.)/4.)+(Z/2.)))**(-1.)

# Stuff
Pressure = []
radius = []
Luminosity = []
Temperature = []
density = []
opacity = []
Epsilon = []
Mass = []
Ratio = []

# Center
# Mass
m = 0
Mass.append(m)
# Radius
r = 0
radius.append(r)
# Luminosity
L = 0
Luminosity.append(L)
# Pressure
Pc = 10**(16.675)
Pressure.append(Pc)
# Temperature
Tc = 10**(7.46)
Temperature.append(Tc)
# Density
rho_c = ((Pc - ((1./3.)*a*(Tc**4.)))*mu)/(NA*k*Tc)
density.append(rho_c)
# Opacity
logR = [-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5]
logT = np.genfromtxt('OPAL_X0.7_Y0.28.dat',skip_header=1,usecols=(0),unpack=True)
logK = np.genfromtxt('OPAL_X0.7_Y0.28.dat', skip_header=1, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14))
f = interpolate.interp2d(logR, logT, logK)
kappa_c = 10**(f(np.log10(rho_c/((Tc/(10**6))**3)), np.log10(Tc)))
opacity.append(kappa_c[0])
# Energy generation
epsilon_pp = ((2.4*10**4)*rho_c*(X**2)*np.exp(-3.38/((Tc/10**9)**(1./3.))))/((Tc/10**9)**(2./3.))	# pp chain
epsilon_CNO = ((4.4*10**25)*rho_c*X*Z*np.exp(-15.228/((Tc/10**9)**(1./3.))))/((Tc/10**9)**(2./3.))	# CNO cycle
epsilon_3alpha = ((5.*10**8)*(rho_c**2)*(Y**3)*np.exp(-4.4/(Tc/10**9)))/((Tc/10**9)**3)	# 3alpha process
epsilon_c = epsilon_pp + epsilon_CNO + epsilon_3alpha
Epsilon.append(epsilon_c)
# Adiabatic index
Pgas = (NA*k*rho_c*Tc)/mu
Beta = Pgas/Pc
Gamma2_c = (32.-(24.*Beta)-(3.*(Beta**2)))/(24.-(18.*Beta)-(3.*Beta**2))


# Offset
# Mass
m = dm
Mass.append(m)
# Radius
r = ((3.*m)/(4.*math.pi*rho_c))**(1./3.)
radius.append(r)
# Luminosity
L = epsilon_c*m
Luminosity.append(L)
# Pressure
P = Pc - ((2./3.)*math.pi*G*(rho_c**2)*(r**2))
Pressure.append(P)
# Temperature
LHS = L
RHS = ((16.*math.pi*a*c*G*(1.-(1./Gamma2_c))*(Tc**4.)*m)/(3.*kappa_c[0]*Pc))
Ratio.append(LHS/RHS)
if LHS > RHS :
	T = Tc - ((1-(1./Gamma2_c))*((2*math.pi*G*(rho_c**2)*Tc*(r**2))/Pc))	# convective
else:
	T = Tc - ((kappa_c*(rho_c**2)*epsilon_c*(r**2))/(8*a*c*(Tc**3)))	# radiative
Temperature.append(T)

# Density
rho = ((P - ((1./3.)*a*(T**4.)))*mu)/(NA*k*T)
density.append(rho)
# Opacity
logR = [-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5]
logT = np.genfromtxt('OPAL_X0.7_Y0.28.dat',skip_header=1,usecols=(0),unpack=True)
logK = np.genfromtxt('OPAL_X0.7_Y0.28.dat', skip_header=1, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14))
f = interpolate.interp2d(logR, logT, logK)
kappa = 10**(f(np.log10(rho/((T/(10**6))**3)), np.log10(T)))
opacity.append(kappa[0])
# Energy generation
epsilon_pp = ((2.4*10**4)*rho*(X**2)*np.exp(-3.38/((T/10**9)**(1./3.))))/((T/10**9)**(2./3.))	# pp chain
epsilon_CNO = ((4.4*10**25)*rho*X*Z*np.exp(-15.228/((T/10**9)**(1./3.))))/((T/10**9)**(2./3.))	# CNO cycle
epsilon_3alpha = ((5.*10**8)*(rho**2)*(Y**3)*np.exp(-4.4/(T/10**9)))/((T/10**9)**3)	# 3alpha process
epsilon = epsilon_pp + epsilon_CNO + epsilon_3alpha
Epsilon.append(epsilon)
# Adiabatic index
Pgas = (NA*k*rho*T)/mu
Beta = Pgas/P
Gamma2 = (32.-(24.*Beta)-(3.*(Beta**2)))/(24.-(18.*Beta)-(3.*Beta**2))

dr = dm/(4.*math.pi*(r**2.)*rho)
dP = (-G*m*dm)/(4*math.pi*(r**4.))
dL = epsilon*dm
LHS = L
RHS = ((16.*math.pi*a*c*G*(1.-(1./Gamma2))*(T**4.)*m)/(3.*kappa[0]*P))
Ratio.append(LHS/RHS)
dT_r = (-3.*kappa[0]*L*dm)/(64.*(math.pi**2.)*a*c*(r**4.)*(T**3.))	# radiation
dT_c = -(1.-(1./Gamma2))*((G*m*T*dm)/(4.*math.pi*(r**4.)*P))	# convection

print(dr, dP, dL, dT_c, LHS/RHS)

# Integrate the equations of stellar structure outwards to surface using microphysics
# while m < M:
# 	dr = dm/(4.*math.pi*(r**2.)*rho)
# 	dP = (-G*m*dm)/(4*math.pi*(r**4.))
# 	dL = epsilon*dm
# 	LHS = L
# 	RHS = ((16.*math.pi*a*c*G*(1.-(1./Gamma2))*(T**4.)*m)/(3.*kappa[0]*P))
# 	Ratio.append(LHS/RHS)
# 	dT_r = (-3.*kappa[0]*L*dm)/(64.*(math.pi**2.)*a*c*(r**4.)*(T**3.))	# radiation
# 	dT_c = -(1.-(1./Gamma2))*((G*m*T*dm)/(4.*math.pi*(r**4.)*P))	# convection
# 	P = P + dP
# 	if P<0:
# 		break
# 	if LHS > RHS :
# 		T = T + dT_c
# 	else:
# 		T = T + dT_r
# 	if T<0:
# 		break
# 	r = r + dr
# 	L = L + dL
# 	m = m + dm
# 	rho = ((P - ((1./3.)*a*(T**4.)))*mu)/(NA*k*T)
# 	Pgas = (NA*k*rho*T)/mu
# 	Beta = Pgas/P
# 	Gamma2 = (32.-(24.*Beta)-(3.*(Beta**2)))/(24.-(18.*Beta)-(3.*Beta**2))
# 	epsilon_pp = ((2.4*10**4)*rho*(X**2)*np.exp(-3.38/((T/10**9)**(1./3.))))/((T/10**9)**(2./3.))	# pp chain
# 	epsilon_CNO = ((4.4*10**25)*rho*X*Z*np.exp(-15.228/((T/10**9)**(1./3.))))/((T/10**9)**(2./3.))	# CNO cycle
# 	epsilon_3alpha = ((5.*10**8)*(rho**2)*(Y**3)*np.exp(-4.4/(T/10**9)))/((T/10**9)**3)	# 3alpha process
# 	epsilon = epsilon_pp + epsilon_CNO + epsilon_3alpha
# 	kappa = 10**(f(np.log10(rho/((T/(10**6))**3)), np.log10(T)))
# 	Mass.append(m)
# 	Pressure.append(P)
# 	Temperature.append(T)
# 	radius.append(r)
# 	Luminosity.append(L)
# 	density.append(rho)
# 	Epsilon.append(epsilon)
# 	opacity.append(kappa[0])

# Mass = np.array(Mass)
# radius = np.array(radius)
# Luminosity = np.array(Luminosity)
# Ratio = np.array(Ratio)
# Pressure = np.array(Pressure)
# Temperature = np.array(Temperature)
# density = np.array(density)
# opacity = np.array(opacity)
# Epsilon = np.array(Epsilon)

# Mass = Mass/Msun


# # Plot log(Pressure) vs m(Msun) (P in cgs)
# fig1 = plt.figure(1)
# plt.plot(Mass,np.log10(Pressure))
# plt.xlabel('Mass in Msun')
# plt.ylabel('log(Pressure)')
# plt.title('log(Pressure) vs. Mass')

# # Plot log(Temperature) vs m(Msun) (T in Kelvin)
# fig2 = plt.figure(2)
# plt.plot(Mass,np.log10(Temperature))
# plt.xlabel('Mass in Msun')
# plt.ylabel('log(Temperature)')
# plt.title('log(Temperature) vs. Mass')

# # Plot log(density) vs m(Msun) (density in cgs)
# fig3 = plt.figure(3)
# plt.plot(Mass,np.log10(density))
# plt.xlabel('Mass in Msun')
# plt.ylabel('log(density)')
# plt.title('log(density) vs. Mass')

# # Plot radius vs m(Msun) (r in solar radii)
# fig4 = plt.figure(4)
# plt.plot(Mass,(radius/Rsun))
# plt.xlabel('Mass in Msun')
# plt.ylabel('r/Rsun')
# plt.title('Radius vs. Mass')

# # Plot luminosity(Lsun) vs m(Msun)
# fig5 = plt.figure(5)
# plt.plot(Mass,(Luminosity/Lsun))
# plt.xlabel('Mass in Msun')
# plt.ylabel('L/Lsun')
# plt.title('Luminosity vs. Mass')

# # Plot opacity vs m(Msun) (in cgs)
# fig6 = plt.figure(6)
# plt.plot(Mass,opacity)
# plt.xlabel('Mass in Msun')
# plt.axis([0.,6.5,0.,10.])
# plt.ylabel('Opacity')
# plt.title('Opacity vs. Mass')

# # Plot energy generation rate vs m(Msun) (cgs)
# fig7 = plt.figure(7)
# plt.plot(Mass,Epsilon)
# plt.xlabel('Mass in Msun')
# plt.ylabel('Epsilon')
# plt.title('Energy Generation Rate vs. Mass')

# # Plot ratio of radiative to adiabatic gradient vs m(Msun) (>1 in convective regions)
# fig8 = plt.figure(8)
# plt.plot(Mass,Ratio)
# plt.plot((0.,7.),(1.,1.), 'r')
# plt.xlabel('Mass in Msun')
# plt.ylabel('Radiative/Adiabatic gradient')
# plt.axis([0.,6.5,0.,10])
# plt.title('Radiative/Adiabatic gradient vs. Mass')


# plt.show()

import cantera as ct
import numpy as np
import matplotlib.pylab as plt

gas: ct.Solution = ct.Solution('gri30.yaml')

NPOINTS: int = 50
TMIN: int = 300
TMAX: int = 3500

data: dict = {
    'reactants': [
        { 'name': 'CH4', 'X': 1/10.52 },
        { 'name': 'O2', 'X': 2/10.52 },
        { 'name': 'N2', 'X': 2*3.76/10.52 }
    ],
    'products': [
        { 'name': 'CO2', 'X': 1/10.52 },
        { 'name': 'H2O', 'X': 2/10.52 },
        { 'name': 'N2', 'X': 7.52/10.52 }
    ],
    'H_GAS_REAC':       np.zeros(NPOINTS),
    'H_GAS_PROD':       np.zeros(NPOINTS),
    'T':                np.linspace(TMIN, TMAX, NPOINTS)
}

for i in range(NPOINTS):

    X: np.array = np.zeros(gas.n_species)
    for reactant in data['reactants']:
        X[gas.species_index(reactant['name'])] = reactant['X']
    gas.TPX = data['T'][i], 101325, X
    data['H_GAS_REAC'][i] = gas.enthalpy_mass

    X: np.array = np.zeros(gas.n_species)
    for product in data['products']:
        X[gas.species_index(product['name'])] = product['X']
    gas.TPX = data['T'][i], 101325, X
    data['H_GAS_PROD'][i] = gas.enthalpy_mass


"""Estimation of the adiabatic flame temperature"""

Ti: float = 298
Hr: float = np.interp(Ti, data['T'], data['H_GAS_REAC'])
Hp: float = Hr
print(f'Enthalpy of the reactants at {Ti} K: {Hr} J/kg')

Tad: float = np.interp(Hp, data['H_GAS_PROD'], data['T'])
print(f'Adiabatic flame temperature: {Tad} K')


"""Plot of the enthalpy of the reactants and products as a function of the temperature"""

fig, ax1 = plt.subplots()

# Plot of the enthalpy of the reactants
ax1.plot(data['T'], data['H_GAS_REAC']/1e+6, 'b-', label='$H_{r}$')
ax1.plot(data['T'], data['H_GAS_PROD']/1e+6, 'r--', label='$H_{p}$')

# Plot of the adiabatic and initial flame temperature
ax1.scatter(Tad, Hp/1e+6, marker='o', color='r', label= '$T_{ad}$ = ' + str(int(Tad)) + ' K')
ax1.scatter(Ti, Hr/1e+6, marker='o', color='b', label= '$T_{i}$ = ' + str(int(Ti)) + ' K')

ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('Enthalpy [MJ/kg]')
ax1.legend()
ax1.grid(True)
plt.xlim((TMIN,TMAX))
plt.show()
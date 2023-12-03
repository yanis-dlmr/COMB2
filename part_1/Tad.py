import cantera as ct
import numpy as np
import matplotlib.pylab as plt

gas: ct.Solution = ct.Solution('gri30.yaml')

NPOINTS = 100
ALPHA_MIN = 0
ALPHA_MAX = 1
PHI_MIN = 0.6
PHI_MAX = 1.0
NPHI = int((PHI_MAX-PHI_MIN)/0.1 + 1)
TI: float = 298

data: dict = {
    'T_ad':     np.array( [np.zeros (NPOINTS) for i in range(NPHI)] ),
    'alpha':    np.linspace(ALPHA_MIN, ALPHA_MAX, NPOINTS),
    'phi':      np.linspace(PHI_MIN, PHI_MAX, NPHI)
}

for i, phi in enumerate(data['phi']):
    for j, alpha in enumerate(data['alpha']):
        X_reac: str = f"CH4:{1-alpha}, H2:{alpha}, O2:{(2-3/2*alpha)/phi}, N2:{3.76*(2-3/2*alpha)/phi}"
        X_prod: str = f"CO2:{1-alpha}, H2O:{2-alpha}, N2:{3.76*(2-3/2*alpha)/phi}, O2:{(2-3/2*alpha)*((1/phi)-1)}"
        
        gas.TPX = TI, 101325, X_reac
        Hr: float = gas.enthalpy_mass

        gas.HPX = Hr, 101325, X_prod
        Tad: float = gas.T

        data['T_ad'][i,j] = Tad


"""Plot of the adiabatic flame temperature as a function of the ratio of H2 to CH4"""

fig, ax1 = plt.subplots()
for i, phi in enumerate(data['phi']):
    ax1.plot(data['alpha'], data['T_ad'][i], label=f'phi = {phi}')
plt.xlabel('alpha [mol$_{H2}$/mol$_{CH4}$]')
plt.ylabel('T$_{ad}$ [K]')
plt.xlim(ALPHA_MIN, ALPHA_MAX)
plt.legend()
plt.show()
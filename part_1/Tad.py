import cantera as ct
import numpy as np
import matplotlib.pylab as plt

from src.graph_tool import Graph_1D

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

# save data
np.save('part_1/data.npy', data)


"""Plot of the adiabatic flame temperature as a function of the ratio of H2 to CH4"""

graph_1d = Graph_1D(fontsize=12, figsize=(6,3))
graph_1d.setup_axis(
    xlabel='$\\alpha$ [-]', 
    ylabel='T$_{ad}$ [K]',
    xmin=ALPHA_MIN,
    xmax=ALPHA_MAX,
    ymin=1600,
    ymax=2600,
    tick=10,
    sci=False
)
colors: list = ['chartjs_red', 'chartjs_orange', 'chartjs_green', 'chartjs_blue', 'chartjs_purple']
for i, phi in enumerate(data['phi']):
    graph_1d.plot(x=data['alpha'], y=data['T_ad'][i], label=f'$\\phi$ = {phi}', marker='', color=colors[i])
graph_1d.save('part_1/T_ad.png', ncol=3, dx=0.05, dy=1.30)

"""Plot phi for a fixed alpha = 0.757"""

NPOINTS = 100
ALPHA = 0.757
PHI_MIN = 0.6
PHI_MAX = 1.0
NPHI = int((PHI_MAX-PHI_MIN)/0.1 + 1)
TI: float = 298

data: dict = {
    'T_ad':     np.zeros(NPHI),
    'phi':      np.linspace(PHI_MIN, PHI_MAX, NPHI)
}

for i, phi in enumerate(data['phi']):
    alpha: float = ALPHA
    X_reac: str = f"CH4:{1-alpha}, H2:{alpha}, O2:{(2-3/2*alpha)/phi}, N2:{3.76*(2-3/2*alpha)/phi}"
    X_prod: str = f"CO2:{1-alpha}, H2O:{2-alpha}, N2:{3.76*(2-3/2*alpha)/phi}, O2:{(2-3/2*alpha)*((1/phi)-1)}"
    
    gas.TPX = TI, 101325, X_reac
    Hr: float = gas.enthalpy_mass

    gas.HPX = Hr, 101325, X_prod
    Tad: float = gas.T

    data['T_ad'][i] = Tad

graph_1d = Graph_1D(fontsize=12, figsize=(5,3))
graph_1d.setup_axis(
    xlabel='$\\phi$ [-]',
    ylabel='T$_{ad}$ [K]',
    xmin=PHI_MIN,
    xmax=PHI_MAX,
    ymin=1600,
    ymax=2600,
    tick=10,
    sci=False
)
graph_1d.plot(x=data['phi'], y=data['T_ad'], marker='', color='chartjs_red')
graph_1d.save(f'part_1/T_ad_alpha_{ALPHA}.png', ncol=3, dy=1)
import cantera as ct
import numpy as np
import matplotlib.pylab as plt

import src.helper.utils as utils

gas: ct.Solution = ct.Solution('gri30.yaml')

NPOINTS: int = 50
TMIN: int = 300
TMAX: int = 3500

species_names: list = ['CH4','O2','CO2','H2O','N2','OH', 'AR']

data: dict = {
    'species_names':    np.array(species_names),
    'CP_list':          np.zeros((len(species_names), NPOINTS)),
    'T':                np.array([ TMIN + (TMAX-TMIN)*i/(NPOINTS-1) for i in range(NPOINTS) ]),
    'indexes':          np.array([ gas.species_index(name) for name in species_names ])
}

for i in range(NPOINTS):
    for j, index in enumerate(data['indexes']):
        X = np.zeros(gas.n_species)
        X[index] = 1
        gas.TPX = data['T'][i], 101325, X
        data['CP_list'][j, i] = gas.cp_mass


"""Plot of the Cp of the mixture as a function of the temperature"""

fig, ax1 = plt.subplots()
for i in range(len(data['indexes'])):
    ax1.plot(data['T'], data['CP_list'][i]/1000, '-', label=utils.species_name_to_label(data['species_names'][i]))
ax1.set_xlabel('$T [K]$')
ax1.set_ylabel('Cp [kJ/kg/K]', color='k')
ax1.legend()
plt.grid()
plt.xlim((TMIN,TMAX))
#plt.show()

"""Generate Air Cp: sum of the Cp of the species weighted by their mole fraction"""

#air_composition: dict = {
#    'N2': 0.7812,
#    'O2': 0.2095,
#    'AR': 0.0093
#}
#
#total_cp: np.array = np.zeros(NPOINTS)
#for species_name, X in air_composition.items():
#    index: int = species_names.index(species_name)
#    total_cp += data['CP_list'][index]*X
#
#fig, ax1 = plt.subplots()
#ax1.plot(data['T'], total_cp/1000, '-', label='Air (sum)')

"""Generate Air Cp with gas.cp_mass"""

#air: ct.Solution = ct.Solution('gri30.yaml')
#
#species: list = ['N2','O2','AR']
#
#
#data:dict = {
#    'species_names':    ['N2','O2','AR'],
#    'T':                np.array([ TMIN + (TMAX-TMIN)*i/(NPOINTS-1) for i in range(NPOINTS) ]),
#    'CP_list':          np.zeros(NPOINTS),
#    'indexes':          np.array([ air.species_index(name) for name in ['N2','O2','AR'] ])
#}
#
#for i in range(NPOINTS):
#    X = np.zeros(air.n_species)
#    X[data['indexes']] = [0.7812, 0.2095, 0.0093]
#    air.TPX = data['T'][i], 101325, X
#    data['CP_list'][i] = air.cp_mass
#    
#ax1.plot(data['T'], data['CP_list']/1000, '-', label='Air (full gas))')

"""Generate Air Cp with gas.cp_mass"""

composition_air: str = "N2:0.7812, O2:0.2095, AR:0.0093"
composition_reac: str = "CH4:1, N2:7.52, O2:2"
air: ct.Solution = ct.Solution('gri30.yaml')
data: dict = {
    'CP_list_air':          np.zeros(NPOINTS),
    'CP_list_reac':         np.zeros(NPOINTS),
    'T':                np.array([ TMIN + (TMAX-TMIN)*i/(NPOINTS-1) for i in range(NPOINTS) ])
}

for i in range(NPOINTS):
    air.TPX = data['T'][i], 101325, composition_air
    data['CP_list_air'][i] = air.cp_mass
    air.TPX = data['T'][i], 101325, composition_reac
    data['CP_list_reac'][i] = air.cp_mass

ax1.plot(data['T'], data['CP_list_air']/1000, '--', label='Air (full gas)')
ax1.plot(data['T'], data['CP_list_reac']/1000, '--', label='Reactifs (full gas)')

ax1.set_xlabel('$T [K]$')
ax1.set_ylabel('Cp [kJ/kg/K]', color='k')
ax1.legend()
plt.xlim((TMIN,TMAX))
plt.show()
import os, sys

"""
Tutorial n�1: Global combustion properties  
Course : Combustion 2 
Energy Engineering Department
INSA Rouen-Normandie
2023-2024 term

Authors: Pradip Xavier 
"""

import cantera as ct
import numpy as np
import matplotlib.pylab as plt

gas = ct.Solution('gri30.yaml')

# We define the "object" gas which represents a phase: a gaseous mixture whose properties are specified in the GRI3.0 mechanism
# This object is a "Solution" type imported via Cantera (prefix ct.) and it is associated to many functions
# Once loaded, it is possible to display all information regarding gas

#print(gas())                          # print the state of gas
print("\nNumber of species : ", gas.n_species)                  # print species number in GRI3.0 mech.
print("\nList of species :\n", gas.species_names)              # print species names in GRI3.0 mech.
name = 'CH4'
ifuel = gas.species_index(name)      # Identify the index associated to the CH4 species in GRI3.0 mech.
print(f"\nIndex of the species {name} : {ifuel}")                          # print the index of species CH4 in GRI3.0 mech.

#help(gas) # allows to access all modules/functions/attributes of the object gas
#print(gas.species('H').thermo.coeffs)
#print (gas.atomic_weights)          #
#print (gas.T)                       # With the IPython terminal, type 'gas.T'


# Possible to modify the thermodynamics state of an object, by fixing 2 intensive variables of the phase
#gas.TP = -- , --          temperature, pressure
#gas.TD = -- , --          temperature, density
#gas.HP = -- , --          specific enthalpy, pressure
#gas.UV = -- , --          specific internal energy, specific volume
#gas.SP = -- , --          specific entropy, pressure
#gas.SV = -- , --          specific entropy, specific volume


gas.TP = 350, 200000       # T=350 K et P = 200000 Pa
#print(gas())

# The mixture composition can be set as follows :

gas.TPX=350,200000,'CH4:1,O2:2,N2:7.52'      # temperature, pressure, composition
print(gas())

# X deals with molar fractions (Y are mass fractions)
# Formating is: <species names>:<mole numbers> (or mass)
# the mole number (or mass) will be automatically converted into molar fractions (or mass)

# The previous method to impose a composition is not flexible.
# An alternative to set the mixture composition:

X = np.zeros(gas.n_species)          # define a vector filled with zeros. The size is the number of species
ich4 = gas.species_index('CH4')      # Identify the index associated to the CH4 species in GRI3.0 mech.
X[ich4] = 1                          # set the mole number to the methane location in vector X
io2 = gas.species_index('O2')        # Repeat for the other species
X[io2] = 2                          
in2 = gas.species_index('N2')      
X[in2] = 7.52                         

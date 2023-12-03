import cantera as ct
import numpy as np
import matplotlib.pylab as plt
from dataclasses import dataclass

from src.graph_tool import Graph_1D
import src.helper.utils as utils

@dataclass
class Combustion:
    """Study of the reaction CH4 + 2O2 = CO2 + 2H2O
    """
    __GAS: ct.Solution = ct.Solution('gri30.yaml')
    __T: int = 300
    __P: int = 101255
    __X_REACTIFS_1 = {
        'CH4': 1,
        'O2': 2
    }
    __X_REACTIFS_2 = {
        'N2': 0.5,
        'O2': 0.5
    }
    __TEMP_MIN: int = 300
    __TEMP_MAX: int = 3500
    __NPOINTS: int = 50
    
    def __post_init__(self):
        self.__data: dict = {
            'temp':       np.linspace(self.__TEMP_MIN, self.__TEMP_MAX, self.__NPOINTS),
            'xeq':        np.zeros((self.__GAS.n_species, self.__NPOINTS)),
            'tf':         np.zeros(self.__NPOINTS)
        }

    def compute_equilibrium_1(self):
        """Compute the equilibrium of the reaction CH4 + 2O2 = CO2 + 2H2O
        """
        print('\nEquilibrium of the reaction CH4 + 2O2 = CO2 + 2H2O')
        self.__GAS.TPX = self.__T, self.__P, self.__X_REACTIFS_1
        self.__GAS.equilibrate('HP')
        print(f'Final temperature = {self.__GAS.T:.1f} K')
        print(f'Molar fraction of CH4 = ', self.__GAS["CH4"].X)

    def compute_equilibrium_2(self):
        """Compute the equilibrium of the reaction 0.5N2 + 0.5O2 = (1-alpha)NO + alphaN2 + 0.5O2
        """
        print('\nEquilibrium of the reaction 0.5N2 + 0.5O2 = (1-alpha)NO + alphaN2 + 0.5O2')
        self.__GAS.TPX = self.__T, self.__P, self.__X_REACTIFS_2
        self.__GAS.equilibrate('HP')
        print(f'Final temperature = {self.__GAS.T:.1f} K')
        print(f'Molar fraction of NO = ', self.__GAS["NO"].X)
        print(f'Alpha = ', 1 - self.__GAS["NO"].X)
        print(f'Molar fraction of N2 = ', self.__GAS["N2"].X)
        print(f'Alpha = ', 2 * self.__GAS["N2"].X)
        print(f'Molar fraction of O2 = ', self.__GAS["O2"].X)
        print(f'Alpha = ', 2 * self.__GAS["O2"].X)

        print(f'Gibbs energy = {self.__GAS.gibbs_mole:.1f} J/mol')

    def compute_evolution(self):
        """Loop on temperature
        """
        for i in range(self.__NPOINTS):
            self.__GAS.TPX = self.__data['temp'][i], self.__P, self.__X_REACTIFS_2
            self.__GAS.equilibrate('HP')
            self.__data['xeq'][:, i] = self.__GAS.X
            self.__data['tf'][i] = self.__GAS.T

    def plot(self):
        """Plot the evolution of the molar fraction of each species
        """
        graph_1d = Graph_1D()
        
        graph_1d.setup_axis(xlabel='Temperature [K]', ylabel='Molar fraction', xmin=self.__TEMP_MIN, xmax=self.__TEMP_MAX, ymin=0, ymax=1)
        for species in self.__GAS.species_names:
            graph_1d.plot(self.__data['temp'], self.__data['xeq'][self.__GAS.species_index(species), :], label=utils.species_name_to_label(species), marker='')
        
        graph_1d.add_axis()
        graph_1d.setup_secondary_axis(ylabel='Final temperature [K]', ymin=300, ymax=3500, color='chartjs_red')
        graph_1d.plot(self.__data['temp'], self.__data['tf'], label='Final temperature', color='chartjs_red', marker='', axis_number=1)
        
        graph_1d.show()
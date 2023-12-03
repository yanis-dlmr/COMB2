import cantera as ct
import numpy as np
import matplotlib.pylab as plt
from dataclasses import dataclass

from src.graph_tool import Graph_1D
import src.helper.utils as utils

@dataclass
class Gibbs_NO:
    __GAS: ct.Solution = ct.Solution('gri30.yaml')
    __P0: int = 101325
    __T2: int = 3000
    __R: float = 8.314
    __NALPHA: int = 100
    __ALPHA_MIN: float = 0.00001
    __ALPHA_MAX: float = 0.99999
    __PROD_SPECIES = ['NO', 'N2', 'O2']
    
    def __post_init__(self):
        self.__data: dict = {
            'alpha':    np.linspace(self.__ALPHA_MIN, self.__ALPHA_MAX, self.__NALPHA),
            'G':    {
                'total':    np.zeros(self.__NALPHA),
                **{species: np.zeros(self.__NALPHA) for species in self.__PROD_SPECIES}
            }
        }
    
    @staticmethod
    def coefficient(alpha: float) -> list:
        """Return the coefficient of the species in the reaction
        """
        return [1 - alpha, 0.5*alpha, 0.5*alpha]

    def compute(self):
        """Compute the data
        """
        for i in range(self.__NALPHA):
            for j, species in enumerate(self.__PROD_SPECIES):
                coefficents: list = self.coefficient(self.__data['alpha'][i])
                Px: float = coefficents[j] * self.__P0
                self.__GAS.TPX = self.__T2, Px, {species: 1}
                self.__data['G'][species][i] = self.__GAS.gibbs_mole
                self.__data['G']['total'][i] += coefficents[j] * (self.__data['G'][species][i] + self.__R * self.__T2 * np.log(Px))

        for i, alpha in enumerate(self.__data['alpha']):
            self.__GAS.TPX = self.__T2, self.__P0, {'NO': 1 - alpha, 'N2': alpha*0.5, 'O2': alpha*0.5}
            self.__data['G']['total'][i] = self.__GAS.gibbs_mole

    def plot(self):
        """Plot the data
        """
        graph_1d = Graph_1D()
        graph_1d.setup_axis(xlabel='Dissociated molar fraction of nitric oxide, alpha', ylabel='G [kJ/kmol]', xmin=self.__ALPHA_MIN, xmax=self.__ALPHA_MAX, ymin=-7.6e5, ymax=-6.7e5, sci=False)
        graph_1d.plot(x=self.__data['alpha'], y=self.__data['G']['total']/1000, color='chartjs_blue', label='Evolution of the Gibbs free energy of the mixture', marker='')
        #for species in self.__PROD_SPECIES:
        #    graph_1d.plot(x=self.__data['alpha'], y=self.__data['G'][species]/1000, color='chartjs_red', label=f'Gibbs free energy of {species}', marker='')
        graph_1d.show(dx=0.25)
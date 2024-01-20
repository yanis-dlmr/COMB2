import numpy as np
import cantera as ct
from dataclasses import dataclass

from src.graph_tool import Graph_1D
import src.helper.utils as utils

R: float = 8.31446261815324

@dataclass()
class Simulation:
    """Class to store the data of a simulation and do the computation
    """
    __gas: ct.Solution = ct.Solution('gri30.yaml')
    __T0: float = 300.0
    __P0: float = 101325.0
    __ALPHA_MIN: float = 0.0
    __ALPHA_MAX: float = 1.0
    __NALPHA: int = 11
    __PHI_MIN: float = 0.6
    __PHI_MAX: float = 1.0
    __NPHI: int = 5
    __data: dict = None

    def __post_init__(self):
        """Initialize the data dictionnary
        """
        self.__data = {
            'alpha':            np.linspace(self.__ALPHA_MIN, self.__ALPHA_MAX, self.__NALPHA),
            'phi':              np.linspace(self.__PHI_MIN, self.__PHI_MAX, self.__NPHI),
            'Tf':               np.zeros((self.__NALPHA, self.__NPHI)),
            'species_names':    self.__gas.species_names,
            'X_PROD':           np.zeros((self.__gas.n_species, self.__NALPHA, self.__NPHI))
        }

    @staticmethod
    def get_X_reactifs(alpha: float, phi: float) -> dict:
        """Return the molar fraction of the reactants for the reaction Methan/Hydrogen/Air
        """
        return {
            "CH4": 1 - alpha,
            "H2": alpha,
            "O2": (2 - 3/2 * alpha) / phi,
            "N2": 3.76 * (2 - 3/2 * alpha) / phi
        }
        
    @staticmethod
    def arrhenius(T: float, A: float, b: float, E: float) -> float:
        """Return the Arrhenius reaction rate
        """
        return A * T**b * np.exp(-E / (T))

    def compute(self):
        """Compute the data
        """
        for i, alpha in enumerate(self.__data['alpha']):
            for j, phi in enumerate(self.__data['phi']):
                X_REACTIFS: dict = self.get_X_reactifs(alpha=alpha, phi=phi)
                self.__gas.TPX = self.__T0, self.__P0, X_REACTIFS
                self.__gas.equilibrate('HP')
                self.__data['Tf'][i, j] = self.__gas.T
                self.__data['X_PROD'][:, i, j] = self.__gas.X

    def get_results(self) -> dict:
        """
        Return the results
        """
        return self.__data
    
    def plot_T_f_1(self):
        """
        Plot T_f = f(alpha) for phi = 1
        """
        graph = Graph_1D()
        graph.setup_axis(xlabel='$\\alpha$ [-]', ylabel='$T_f$ [K]', xmin=self.__ALPHA_MIN, xmax=self.__ALPHA_MAX, ymin=2200, ymax=2400, tick=10, sci=False)
        colors: list = ['chartjs_red', 'chartjs_orange', 'chartjs_green', 'chartjs_blue', 'chartjs_purple']
        for i, phi in enumerate(self.__data['phi']):
            if phi == 1:
                graph.plot(x=self.__data['alpha'], y=self.__data['Tf'][:, i], label=f'$\phi = {phi}$', color=colors[i])
        graph.save(filename='part_2/results/T_f_1', dx=0.4, dy=1.12)

    def plot_T_f(self):
        """
        Plot T_f = f(alpha)
        """
        graph = Graph_1D()
        graph.setup_axis(xlabel='$\\alpha$ [-]', ylabel='$T_f$ [K]', xmin=self.__ALPHA_MIN, xmax=self.__ALPHA_MAX, ymin=1600, ymax=2400, tick=10, sci=False)
        colors: list = ['chartjs_red', 'chartjs_orange', 'chartjs_green', 'chartjs_blue', 'chartjs_purple']
        for i, phi in enumerate(self.__data['phi']):
            graph.plot(
                x=self.__data['alpha'], y=self.__data['Tf'][:, i], label=f'$\phi = {phi}$', color=colors[i])
        graph.save(filename='part_2/results/T_f', ncol=3, dx=0.19, dy=1.17)
    
    def compare_T_f(self):
        """
        Compare T_f = f(alpha) with T_ad = f(alpha) for different phi
        """
        data = np.load('part_1/data.npy', allow_pickle=True).item()
        graph = Graph_1D(fontsize=12, figsize=(11, 5))
        graph.setup_axis(xlabel='$\\alpha$ [-]', ylabel='T [K]', xmin=self.__ALPHA_MIN, xmax=self.__ALPHA_MAX, ymin=1600, ymax=2600, sci=False)
        colors: list = ['chartjs_red', 'chartjs_orange', 'chartjs_green', 'chartjs_blue', 'chartjs_purple']
        for i, phi in enumerate(self.__data['phi']):
            graph.plot(x=self.__data['alpha'], y=self.__data['Tf'][:, i], label=f'$\\phi$ = {phi} T$_f$', marker='o', color=colors[i], linestyle='-.')
            graph.plot(x=data['alpha'], y=data['T_ad'][i], label=f'$\\phi$ = {phi} ' + 'T$_{ad}$', marker='', color=colors[i], linestyle=':')
        graph.save(filename='part_2/results/T_f_compare', ncol=5, dx=0, dy=1.19)
        

    def plot_X_PROD(self):
        """
        Plot X_PROD = f(alpha) for different phi
        """
        #for i, phi in enumerate(self.__data['phi']):
        #    graph = Graph_1D()
        #    graph.setup_axis(
        #        xlabel='$\\alpha$ [-]', ylabel='$X_{PROD}$ [-]', xmin=self.__ALPHA_MIN, xmax=self.__ALPHA_MAX, ymin=0, ymax=1, tick=10, sci=False)
        #    for j, species in enumerate(self.__data['species_names']):
        #        if (np.max(self.__data['X_PROD'][j, :, i]) > 1e-3):
        #            graph.plot(x=self.__data['alpha'], y=self.__data['X_PROD'][j, :, i], label=f'{utils.species_name_to_label(species)}', marker='')
        #    graph.show(dx=0.25, ncol=4)
        
        colors: list = ['chartjs_red', 'chartjs_orange', 'chartjs_green', 'chartjs_blue', 'chartjs_purple']
        for i, species in enumerate(self.__data['species_names']):
            if (np.max(self.__data['X_PROD'][i, :, :]) > 1e-6):
                graph = Graph_1D(fontsize=12, figsize=(6, 4))
                graph.setup_axis(
                    xlabel='$\\alpha$ [-]', ylabel=f'{utils.species_name_to_label(species)} [ppm]', 
                    xmin=self.__ALPHA_MIN, 
                    xmax=self.__ALPHA_MAX, 
                    #ymin=np.min(self.__data['X_PROD'][i, :, :]) * 1e6, 
                    #ymax=np.max(self.__data['X_PROD'][i, :, :]) * 1e6,
                    sci=False
                )
                for j, phi in enumerate(self.__data['phi']):
                    graph.plot(x=self.__data['alpha'], y=self.__data['X_PROD'][i, :, j] * 1e6, label=f'$\phi = {phi}$', marker='o', color=colors[j])
                #graph.show(dx=0.25, ncol=3)
                graph.save(filename='part_2/species/' + species + '.png', ncol=3, dx=0.05, dy=1.2)

    def plot_reaction_rate_constant_N1_N2(self):
        """
        Plot reaction rate constant for reaction (N1)
        """
        T: list = np.array([ T for T in range(1000, 4500) ])
        k_N1_f: list = np.array([ self.arrhenius(_T, 1.8e+11, 0, 38370) for _T in T ])
        k_N2_f: list = np.array([ self.arrhenius(_T, 1.8e+10, 0, 4680) for _T in T ])
        
        graph = Graph_1D(fontsize=12, figsize=(8, 4))
        graph.setup_axis(xlabel='T$_f$ [K]', ylabel='$k_{N,f}$ [m$^3$/mol/s]', xmin=1000, xmax=4500, ymin=0, ymax=1e8, color='chartjs_red')
        graph.plot(x=T, y=k_N1_f, label='$k_{N1}$', marker='', color='chartjs_red')
        graph.add_axis()
        graph.setup_secondary_axis(ylabel='$k_{N,f}$ [m$^3$/mol/s]', axis_number=1, ymin=0, ymax=1e10, color='chartjs_blue')
        graph.plot(x=T, y=k_N2_f, label='$k_{N2}$', marker='', color='chartjs_blue', axis_number=1)
        graph.save(filename='part_2/results/k_N.png', ncol=2, dx=0.4, dy=1.2)
        
        graph = Graph_1D(fontsize=12, figsize=(6, 4))
        graph.setup_axis(xlabel='1000/T$_f$ [K$^{-1}$]', ylabel='$k_{N,f}$ [m$^3$/mol/s]', axis_number=0, log_scale_y=True, xmax=1, xmin=0.3, ymin=1e-6, ymax=1e7)
        graph.plot(x=1000/T, y=k_N1_f, marker='', color='chartjs_red')
        graph.secondary_xaxis(functions=(self.x2T, self.T2x), xlabel='T$_f$ [K]')
        graph.save(filename='part_2/results/k_N1.png', ncol=3, dx=0.05, dy=0.8)
        
        graph = Graph_1D(fontsize=12, figsize=(6, 4))
        graph.setup_axis(xlabel='1000/T$_f$ [K$^{-1}$]', ylabel='$k_{N,f}$ [m$^3$/mol/s]', axis_number=0, log_scale_y=True, xmax=1, xmin=0.3, ymin=1e7, ymax=1e11)
        graph.plot(x=1000/T, y=k_N2_f, marker='', color='chartjs_blue')
        graph.secondary_xaxis(functions=(self.x2T, self.T2x), xlabel='T$_f$ [K]')
        graph.save(filename='part_2/results/k_N2.png', ncol=3, dx=0.05, dy=0.8)
    
    def x2T(self, x: float) -> float:
        """
        Return the temperature of the gas for a given x
        """
        return 1000/x
    
    def T2x(self, T: float) -> float:
        """
        Return the x of the gas for a given temperature
        """
        return 1000/T
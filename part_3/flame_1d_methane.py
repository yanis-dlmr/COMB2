import cantera as ct
import numpy as np
import matplotlib.pylab as plt
from dataclasses import dataclass

from src.graph_tool import Graph_1D

@dataclass
class Simulation:

    __loglevel  = 1
    __refine_grid = True
    __GAS: ct.Solution = ct.Solution('gri30.yaml')
    __T_INIT: int = 293
    __P_INIT: int = 101325
    __PHI: int = 1
    # Grid specfications with refinement at inlet and outlet, 6 points in x-direction
    __GRID_INIT = 5*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'f')/3.0

    def initialize(self, alpha: float = 1.0) -> None:
        """Initialize the data dictionnary
        """
        self.current_alpha = alpha
        self.__GAS.TPX = self.__T_INIT, self.__P_INIT, self.get_X_reactifs(alpha)
        self.__F: ct.FreeFlame = ct.FreeFlame(self.__GAS, self.__GRID_INIT)
        self.__F.flame.set_steady_tolerances(default=[1.0e-5, 1.0e-8])
        self.__F.flame.set_transient_tolerances(default=[1.0e-5, 1.0e-8])
        self.__F.inlet.X = self.get_X_reactifs(alpha)
        self.__F.inlet.T = self.__T_INIT

    def get_X_reactifs(self, alpha: float = 1.0) -> dict:
        """Return the molar fraction of the reactants for the reaction Methan/Air
        """
        phi = 1.0
        return {
            "CH4": 1 - alpha,
            "H2": alpha,
            "O2": (2 - 3/2 * alpha) / phi,
            "N2": 3.76 * (2 - 3/2 * alpha) / phi
        }  

    def compute(self) -> None:
        """Compute the data
        """
        _ratio, _slope, _curve, _prune = 5.0, 0.02, 0.02, 0.02
        self.__F.energy_enabled = True
        self.__F.set_max_jac_age(50, 50)
        self.__F.set_time_step(0.1e-06, [2, 5,10, 20, 80])
        self.__F.set_refine_criteria(ratio = _ratio, slope = _slope, curve = _curve, prune = _prune)
        self.__F.solve(self.__loglevel, self.__refine_grid)

    def compute_flame_speed(self, species: str) -> float:
        """Flame speed s0L as a function of the integral of the reaction rate of the fuel.
        This approach can work with any majority species, present in fresh gases or in burned gases
        If species is freshgases, the flame speed is the speed of the fresh gases.
        """
        if species == "freshgases":
            return self.__F.velocity[0]
        else:
            _reaction_rate: list = self.__F.net_production_rates[self.__GAS.species_index(species)]
            _Y_f_1: float = self.__F.Y[self.__GAS.species_index(species), 0]
            _Y_f_2: float = self.__F.Y[self.__GAS.species_index(species), -1]
            _volumic_mass: float = self.__F.density[0]
            _gas = ct.Solution('gri30.yaml')
            _gas.TPX = self.__T_INIT, self.__P_INIT, {species: 1.0}
            _molar_mass: float = _gas.mean_molecular_weight
            S0l = np.trapz(_reaction_rate, self.__F.grid) / (_Y_f_2 - _Y_f_1) / _volumic_mass * _molar_mass
            return S0l
    
    def compute_old(self) -> None:
        """Compute the data
        """
        _ratio: list = [7.0, 5.0, 5.0, 5.0, 5.0, 5.0]
        _slope: list = [1, 0.5, 0.3, 0.1, 0.05, 0.02]
        _curve: list = _slope
        _prune: list = [0, 0, 0, 0, 0.01, 0.02]
        
        for i, (r, s, c, p) in enumerate(zip(_ratio, _slope, _curve, _prune)):
            if i == 0:
                self.__F.energy_enabled = False
                self.__F.set_max_jac_age(50, 50)
                self.__F.set_time_step(0.1e-06, [2, 5,10, 20, 80])
                self.__F.set_refine_criteria(ratio = r, slope = s, curve = c, prune = p)
                self.__F.solve(self.__loglevel, self.__refine_grid)
                self.__F.energy_enabled = True
            else:
                self.__F.set_refine_criteria(ratio = r, slope = s, curve = c, prune = p)
                self.__F.solve(self.__loglevel, self.__refine_grid)
    
    def get_flame_speed(self) -> float:
        """Return the flame speed
        """
        return self.__F.velocity

    def get_temperature(self) -> np.array:
        """Return the temperature
        """
        return self.__F.T
    
    def get_grid(self) -> np.array:
        """Return the grid
        """
        return self.__F.grid

    def get_net_production_rates_hydrogen(self) -> np.array:
        """Return the net production rates of hydrogen
        """ 
        return self.__F.net_production_rates[self.__GAS.species_index('H')]

    def save_solution(self) -> None:
        """Save the solution
        """
        self.__F.save(f'part_3/results/methane_flame_{self.current_alpha}.yaml', 'solution', 'flame', 'grid', overwrite=True)
    
    def load_solution(self) -> None:
        """Load the solution
        """
        self.__F.restore(f'part_3/results/methane_flame_{self.current_alpha}.yaml', 'solution', 'flame')
    
    def plot_mole_fractions(self) -> None:
        """Plot the mole fractions
        """
        z: np.ndarray = self.get_grid()
        T: np.ndarray = self.get_temperature()
        
        graph_1d = Graph_1D()
        
        graph_1d.setup_axis(xlabel='Axial position [m]', ylabel='Mole fractions [-]', xmin=0.0, xmax=0.05, ymin=0.0, ymax=6e-6, color='chartjs_blue')
        graph_1d.plot(x=z, y=self.__F.X[self.__F.flame.component_index('CH4')], label='$X_{CH4}$', color='chartjs_blue', marker='')
        graph_1d.plot(x=z, y=self.__F.X[self.__F.flame.component_index('O2')], label='$X_{O2}$', color='chartjs_green', marker='' )
        
        graph_1d.add_axis()
        graph_1d.setup_secondary_axis(ylabel='Temperature [K]', ymin=200.0, ymax=2400.0, color='chartjs_red', sci=False, axis_number=1 )
        graph_1d.plot(x=z, y=T, label='Temperature', color='chartjs_red', axis_number=1, linestyle='-', marker='')
        
        graph_1d.show(dx=0.35)
    
    def plot_flow_velocity(self) -> None:
        """Plot the flow velocity
        """
        z: np.ndarray = self.get_grid()
        u: np.ndarray = self.get_flame_speed()
        T: np.ndarray = self.get_temperature()
        
        graph_1d = Graph_1D()
    
        graph_1d.setup_axis(xlabel='Axial position [m]', ylabel='Flow velocity [m/s]', xmin=0.0, xmax=0.05, ymin=0.0, ymax=3, color='chartjs_green')
        graph_1d.plot(x=z, y=u, label='Flow velocity', color='chartjs_green', marker='')
        
        graph_1d.add_axis()
        graph_1d.setup_secondary_axis(ylabel='Temperature [K]', ymin=200.0, ymax=2400.0, color='chartjs_red', sci=False, axis_number=1 )
        graph_1d.plot(x=z, y=T, label='Temperature', color='chartjs_red', marker='', axis_number=1)

        graph_1d.show(dx=0.35)
        
    def plot_flow_density(self) -> None: 
        """Plot the flow density
        """
        z: np.ndarray = self.get_grid()
        T: np.ndarray = self.get_temperature()
        density_mass: np.ndarray = self.__F.density
        
        graph_1d = Graph_1D()
    
        graph_1d.setup_axis(xlabel='Axial position [m]', ylabel='Flow density [kg/m3]', xmin=0.0, xmax=0.05, ymin=0.0, ymax=1.2, color='chartjs_purple')
        graph_1d.plot(x=z, y=density_mass, label='Flow density', color='chartjs_purple', marker='')
    
        graph_1d.add_axis()
        graph_1d.setup_secondary_axis(ylabel='Temperature [K]', ymin=200.0, ymax=2400.0, color='chartjs_red', sci=False, axis_number=1 )
        graph_1d.plot(x=z, y=T, label='Temperature', color='chartjs_red', marker='', axis_number=1)
    
        graph_1d.show(dx=0.35)
    
    def plot_T(self) -> None:
        """Plot the temperature
        """
        z: np.ndarray = self.get_grid()
        T: np.ndarray = self.get_temperature()
        
        graph = Graph_1D(fontsize=12, figsize=(6, 4), cmap='coolwarm')
    
        graph.setup_axis(xlabel='Axial position [m]', ylabel='Temperature [K]', xmin=0.0, xmax=0.05, ymin=200.0, ymax=4000.0, color='chartjs_red')
        graph.plot(x=z, y=T, label='Temperature', color='chartjs_red', marker='')
    
        graph.save(filename='part_3/results/temperature.png', ncol=3, dx=0.23, dy=1.25)
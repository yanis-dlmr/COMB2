import cantera as ct
import numpy as np
import matplotlib.pylab as plt
from dataclasses import dataclass

from src.graph_tool import Graph_1D


@dataclass
class Simulation:

    __loglevel  = 1
    __refine_grid = True
    __GAS: ct.Solution = ct.Solution('gri30.yaml')#, 'gri30_mix')
    __T_INIT: int = 293
    __P_INIT: int = 101325
    __PHI: int = 1
    # Grid specfications with refinement at inlet and outlet, 6 points in x-direction
    __GRID_INIT = 5*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'f')/3.0

    def __post_init__(self):
        """Initialize the data dictionnary
        """
        self.__GAS.TPX = self.__T_INIT, self.__P_INIT, self.get_X_reactifs()
        self.__F: ct.FreeFlame = ct.FreeFlame(self.__GAS, self.__GRID_INIT)
        self.__F.flame.set_steady_tolerances(default=[1.0e-5, 1.0e-8])
        self.__F.flame.set_transient_tolerances(default=[1.0e-5, 1.0e-8])
        self.__F.inlet.X = self.get_X_reactifs()
        self.__F.inlet.T = self.__T_INIT

    def get_X_reactifs(self) -> dict:
        """Return the molar fraction of the reactants for the reaction Methan/Air
        """
        return { "CH4": self.__PHI, "O2": 2, "N2": 3.76 * 2 }

    def compute(self) -> None:
        """Compute the data
        """
        _ratio, _slope, _curve, _prune = 5.0, 0.02, 0.02, 0.02
        self.__F.energy_enabled = True
        self.__F.set_max_jac_age(50, 50)
        self.__F.set_time_step(0.1e-06, [2, 5,10, 20, 80])
        self.__F.set_refine_criteria(ratio = _ratio, slope = _slope, curve = _curve, prune = _prune)
        self.__F.solve(self.__loglevel, self.__refine_grid)
    
    def compute_flame_speed(self, species: str) -> None:
        """Flame speed s0L as a function of the integral of the reaction rate of the fuel.
        This approach can work with any majority species, present in fresh gases or in burned gases
        """
        _reaction_rate: list = self.__F.net_production_rates[self.__GAS.species_index(species)]
        _Y_f_1: float = self.__F.Y[0, self.__GAS.species_index(species)]
        _Y_f_2: float = self.__F.Y[-1, self.__GAS.species_index(species)]
        _volumic_mass: float = self.__F.density[0]
        S0l = np.trapz(_reaction_rate, self.__F.grid) / (_Y_f_2 - _Y_f_1) / _volumic_mass
        print(f"Flame speed for {species} s0L = {S0l} m/s")
        
    
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
    
    def gert_grid(self) -> np.array:
        """Return the grid
        """
        return self.__F.grid
    
    def save_solution(self) -> None:
        """Save the solution
        """
        self.__F.save('methane_flame.yaml', 'solution', 'flame', 'grid', overwrite=True)
    
    def load_solution(self) -> None:
        """Load the solution
        """
        self.__F.restore('part_3/methane_flame.yaml', 'solution', 'flame')
    
    def plot_mole_fractions(self) -> None:
        """Plot the mole fractions
        """
        z: np.ndarray = self.gert_grid()
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
        z: np.ndarray = self.gert_grid()
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
        z: np.ndarray = self.gert_grid()
        T: np.ndarray = self.get_temperature()
        density_mass: np.ndarray = self.__F.density
        
        graph_1d = Graph_1D()
    
        graph_1d.setup_axis(xlabel='Axial position [m]', ylabel='Flow density [kg/m3]', xmin=0.0, xmax=0.05, ymin=0.0, ymax=1.2, color='chartjs_purple')
        graph_1d.plot(x=z, y=density_mass, label='Flow density', color='chartjs_purple', marker='')
    
        graph_1d.add_axis()
        graph_1d.setup_secondary_axis(ylabel='Temperature [K]', ymin=200.0, ymax=2400.0, color='chartjs_red', sci=False, axis_number=1 )
        graph_1d.plot(x=z, y=T, label='Temperature', color='chartjs_red', marker='', axis_number=1)
    
        graph_1d.show(dx=0.35)
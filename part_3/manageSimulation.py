import cantera as ct
import numpy as np
import matplotlib.pylab as plt
from dataclasses import dataclass

from src.graph_tool import Graph_1D
from .flame_1d_methane import Simulation

@dataclass
class ManageSimulation:
    COMPUTE: bool = False
    
    def __post_init__(self) -> None:    
        sim = Simulation()
        self.alphas: np.ndarray = np.linspace(0.0, 1.0, 11)
        
        self.Sl0: dict = {
            "CH4": np.zeros_like(self.alphas),
            "H2": np.zeros_like(self.alphas),
            "N2": np.zeros_like(self.alphas),
            "H2O": np.zeros_like(self.alphas),
            "CO2": np.zeros_like(self.alphas),# it sucks at alpha = 1.0
            "O2": np.zeros_like(self.alphas),
            "CO": np.zeros_like(self.alphas),# it sucks at alpha = 1.0
            "OH": np.zeros_like(self.alphas),
            "NO": np.zeros_like(self.alphas),
            "freshgases": np.zeros_like(self.alphas)
        }
        
        self.temperature: list = []
        self.z: list = []
        self.t_alpha: list = []
        self.net_production_rates_hydrogen: list = []
        
        for _i, _a in enumerate(self.alphas):
            sim.initialize(alpha=_a)
            if self.COMPUTE:
                sim.compute_old()
                sim.save_solution()
            else:
                sim.load_solution()

            for _s in self.Sl0.keys():
                self.Sl0[_s][_i] = sim.compute_flame_speed(_s)
            
            _a = round(_a, 1)
            print(_a, round(_a % 0.2, 1), round(_a % 0.2, 1)%0.2)
            # si i == 0 ou 0.2 ou  0.4 etc tec
            if (_a == 0 or round(_a % 0.2, 1)%0.2 == 0 or _a == 1.0):
                self.temperature.append(sim.get_temperature())
                self.z.append(sim.get_grid())
                self.t_alpha.append(_a)
                
                self.net_production_rates_hydrogen.append(sim.get_net_production_rates_hydrogen())
            
            
            #sim.plot_mole_fractions()
            #sim.plot_flow_velocity()
            #sim.plot_flow_density()
        
        self.plot_T()
        self.plot_net_production_rates_hydrogen()
            
        self.plot_Sl0()
        self.plot_reaction_rate_constant_R3_R84()
    
    def plot_T(self) -> None:
        """Plot the temperature along the flame grid for different alpha coefficients
        """
        colors: list = ['chartjs_red', 'chartjs_orange', 'chartjs_green', 'chartjs_blue', 'chartjs_purple', 'chartjs_pink']
        
        graph_1d = Graph_1D(fontsize=12, figsize=(6, 4))
        graph_1d.setup_axis(xlabel='z [m]', ylabel='T [K]', xmin=0.0, xmax=0.05, sci=False)
        e_list: list = []
        for _i, (_a, _T, _z) in enumerate(zip(self.t_alpha, self.temperature, self.z)):
            graph_1d.plot(x=_z, y=_T, label=f'$\\alpha$ = {_a:.1f}', marker='', color=colors[_i])
            T_gradient: np.ndarray = np.gradient(_T, _z)
            e = (_T[-1] - _T[0]) / np.max(T_gradient) * 1000
            e_list.append(e)
        graph_1d.save(filename='part_3/results/temperature.png', ncol=3, dx=0.1, dy=1.2)
        
        graph_1d_thickness = Graph_1D(fontsize=12, figsize=(6, 4))
        graph_1d_thickness.setup_axis(xlabel='$\\alpha$ [-]', ylabel='e [mm]', xmin=0.0, xmax=1, sci=False)
        graph_1d_thickness.plot(x=self.t_alpha, y=e_list, label=f'', marker='o', color='chartjs_green')
        graph_1d_thickness.save(filename='part_3/results/flame_thickness.png', ncol=3, dx=0.1, dy=1.2)
    
    def plot_net_production_rates_hydrogen(self) -> None:
        """Plot the net production rates of hydrogen along the flame grid for different alpha coefficients
        """
        graph_1d = Graph_1D(fontsize=12, figsize=(6, 4))
        
        graph_1d.setup_axis(xlabel='z [m]', ylabel='$\\dot\\omega_H$ [kmol/m$^3$/s]', xmin=0.018, xmax=0.020, sci=False)
        
        colors: list = ['chartjs_red', 'chartjs_orange', 'chartjs_green', 'chartjs_blue', 'chartjs_purple', 'chartjs_pink']
        
        for _i, (_a, _npr_h, _z) in enumerate(zip(self.t_alpha, self.net_production_rates_hydrogen, self.z)):
            graph_1d.plot(x=_z, y=_npr_h, label=f'$\\alpha$ = {_a:.1f}', marker='', color=colors[_i])
        
        graph_1d.save(filename='part_3/results/net_production_rates_hydrogen.png', ncol=3, dx=0.1, dy=1.2)
            
    def plot_Sl0(self) -> None:
        """Plot the flame speed as a function of the equivalence ratio
        """
        graph_1d = Graph_1D(fontsize=12, figsize=(7, 5))
        
        graph_1d.setup_axis(xlabel='$\\alpha$ [-]', ylabel='Flame speed [cm/s]', xmin=0.0, xmax=1.0, sci=False)
        
        for _s in self.Sl0.keys():
            if _s in ["CO", "CO2"]: #remove last point
                graph_1d.plot(x=self.alphas[:-1], y=self.Sl0[_s][:-1]*100, label=_s)
            else:
                graph_1d.plot(x=self.alphas, y=self.Sl0[_s]*100, label=_s)
        
        graph_1d.save(filename='part_3/results/speed.png', ncol=4, dx=0.0, dy=1.2)
        
        graph_1d = Graph_1D(fontsize=12, figsize=(7, 5))
        graph_1d.setup_axis(xlabel='$\\alpha$ [-]', ylabel='$\\Delta$ Flame speed [cm/s]', xmin=0.0, xmax=1.0, sci=False)
        for _s in self.Sl0.keys():
            reference = self.Sl0["freshgases"]
            if _s != "freshgases":
                if _s in ["CO", "CO2"]: #remove last point
                    graph_1d.plot(x=self.alphas[:-1], y=(self.Sl0[_s][:-1]-reference[:-1])*100, label=_s)
                else:
                    graph_1d.plot(x=self.alphas, y=(self.Sl0[_s]-reference)*100, label=_s)
        graph_1d.save(filename='part_3/results/speed_compare.png', ncol=4, dx=0.1, dy=1.2)

        
    @staticmethod
    def arrhenius(T: float, A: float, b: float, E: float) -> float:
        """Return the Arrhenius reaction rate
        """
        return A * T**b * np.exp(-E / (T * 8.314))

    def plot_reaction_rate_constant_R3_R84(self):
        """
        Plot reaction rate constant for reaction (R3) and (R84)
        """
        T: list = np.array([ T for T in range(1000, 4500) ])
        k_R3_f: list = np.array([ self.arrhenius(_T, 3.87000E+04, 2.7, 6260) for _T in T ])
        k_R84_f: list = np.array([ self.arrhenius(_T, 2.16000E+08, 1.51, 3430) for _T in T ])
        
        graph = Graph_1D(fontsize=12, figsize=(6, 4))
        graph.setup_axis(xlabel='1000/T$_f$ [K$^{-1}$]', ylabel='$k_{N,f}$ [m$^3$/mol/s]', axis_number=0, log_scale_y=True, xmax=1, xmin=0.3)#, ymin=1e-6, ymax=1e7)
        graph.plot(x=1000/T, y=k_R3_f, marker='', color='chartjs_red', label='R3')
        graph.plot(x=1000/T, y=k_R84_f, marker='', color='chartjs_blue', label='R84')
        graph.secondary_xaxis(functions=(self.x2T, self.T2x), xlabel='T$_f$ [K]')
        graph.save(filename='part_3/results/k_R3_R84.png', ncol=3, dx=0.23, dy=1.25)
    
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
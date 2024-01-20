from dataclasses import dataclass
import numpy as np
import cantera as ct

from src.graph_tool import Graph_1D

TIME_MIN: float = 0.0
TIME_MAX: float = 0.0001
NPOINTS: int = 100000
FUEL_FORMATION_ENTHALPY: float = -84e3
O2_FORMATION_ENTHALPY: float = 0
CO2_FORMATION_ENTHALPY: float = -393.5e3
H2O_FORMATION_ENTHALPY: float = -241.8e3
GAS_CONSTANT = 8.31446261815324

@dataclass
class Evolution:
    T0: float = 753
    P0: float = 25.12*ct.one_atm

    time: np.ndarray = np.linspace(TIME_MIN, TIME_MAX, 100)
    temperature: np.ndarray = np.zeros_like(time)
    pressure: np.ndarray = np.zeros_like(time)
    
    concentrations = {
        'F':        np.zeros_like(time),
        'O2':       np.zeros_like(time),
        'PRODUCTS': np.zeros_like(time),
    }
    
    gas: ct.Solution = ct.Solution('gri30.yaml')
    
    def __post_init__(self):
        self.gas.basis = 'molar'
        self.temperature[0] = self.T0
        self.pressure[0] = self.P0
        self.concentrations['F'][0] = 23.049
        self.concentrations['O2'][0] = 80.673
        self.concentrations['PRODUCTS'][0] = 0.0
        
    def get_molar_heat_capacity_fuel(self, T:float, P: float) -> float:
        """Compute the molar heat capacity of the fuel at constant pressure.
        """
        self.gas.TPX = T, P, { 'C2H6': 1 }
        return self.gas.cp / 1000

    def compute(self) -> None:
        """Compute the evolution along the time.
        """
        for _i, _t in enumerate(self.time):
            
            if _i == 0:
                continue
            
            _delta_t: float = _t - self.time[_i-1]
            _prev_temperature: float = self.temperature[_i-1]
            _prev_pressure: float = self.pressure[_i-1]
            
            _dconcentration_fuel_dt : float = - 1.1e12 * np.exp(-15098/_prev_temperature) \
                * (self.concentrations['F'][_i-1])**0.1 * (self.concentrations['O2'][_i-1])**1.65
            self.concentrations['F'][_i] = self.concentrations['F'][_i-1] + _dconcentration_fuel_dt * _delta_t
            
            _molar_reaction_rate_fuel: float = - _dconcentration_fuel_dt
            
            _dconcentration_o2_dt : float = _dconcentration_fuel_dt * (16.66/4.76)
            self.concentrations['O2'][_i] = self.concentrations['O2'][_i-1] + _dconcentration_o2_dt * _delta_t        
    
            _dconcentration_products_dt : float = _molar_reaction_rate_fuel * (5+3.76*16.66/4.76)
            self.concentrations['PRODUCTS'][_i] = self.concentrations['PRODUCTS'][_i-1] + _dconcentration_products_dt * _delta_t
    
            _sum_hformation_reaction_rate: float = FUEL_FORMATION_ENTHALPY * _molar_reaction_rate_fuel
            
            _dT_dt = - _sum_hformation_reaction_rate / (self.get_molar_heat_capacity_fuel(
                T=_prev_temperature,P=_prev_pressure) - GAS_CONSTANT) * _prev_temperature * GAS_CONSTANT / _prev_pressure
            self.temperature[_i] = _prev_temperature + _dT_dt * _delta_t
            
            _dP_dt = _prev_pressure / _prev_temperature * _dT_dt
            self.pressure[_i] = _prev_pressure + _dP_dt * _delta_t

    def plot_concentrations(self) -> None:
        """Plot the concentrations evolution.
        """
        graph_1d = Graph_1D()
        graph_1d.setup_axis(
            xlabel='Time [s]',
            ylabel='Concentration [mol/m3]',
            xmin=TIME_MIN,
            xmax=TIME_MAX,
            ymin=0,
            ymax=100,
            sci=False
        )
        graph_1d.plot(
            x=self.time,
            y=self.concentrations['F'],
            label='Concentration Fuel $C_2H_6$',
            color='chartjs_purple',
            marker='',
        )
        graph_1d.plot(
            x=self.time,
            y=self.concentrations['O2'],
            label='Concentration $O_2$',
            color='chartjs_blue',
            marker='',
        )
        #graph_1d.plot(
        #    x=self.time,
        #    y=self.concentrations['PRODUCTS'],
        #    label='Concentration $PRODUCTS$',
        #    color='chartjs_green',
        #    marker='',
        #)
        graph_1d.add_axis()
        graph_1d.setup_axis(
            xlabel='Time [s]',
            ylabel='Temperature [K]',
            xmin=TIME_MIN,
            xmax=TIME_MAX,
            ymin=0,
            ymax=1000,
            sci=False,
            axis_number=1,
            color='chartjs_red'
        )
        graph_1d.plot(
            x=self.time,
            y=self.temperature,
            label='Temperature',
            color='chartjs_red',
            marker='',
            axis_number=1
        )
        graph_1d.show(dx=0.2)
    
    def plot_pressure(self) -> None:
        """Plot the pressure evolution.
        """
        graph_1d = Graph_1D()
        graph_1d.setup_axis(
            xlabel='Time [s]',
            ylabel='Pressure [atm]',
            xmin=TIME_MIN,
            xmax=TIME_MAX,
            ymin=0,
            ymax=30,
            sci=False
        )
        graph_1d.plot(
            x=self.time,
            y=self.pressure / ct.one_atm,
            label='Pressure',
            color='chartjs_red',
            marker='',
        )
        graph_1d.show(dx=0.2)
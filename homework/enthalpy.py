from dataclasses import dataclass
import numpy as np
import cantera as ct

from src.graph_tool import Graph_1D

TMIN: float = 1
TMAX: float = 3500.0
NPOINTS: int = 100
FUEL_FORMATION_ENTHALPY: float = -84e3

@dataclass
class Enthalpy:
    gas = ct.Solution('gri30.yaml')
    temperature: list = np.linspace(TMIN, TMAX, NPOINTS)
    enthalpy_fuel: list = np.zeros_like(temperature)
    enthalpy_oxidier: list = np.zeros_like(temperature)
    
    def get_molar_heat_capacity_fuel(self, species: dict, T:float, P: float = ct.one_atm) -> float:
        """Compute the molar heat capacity of the fuel at constant pressure.
        """
        self.gas.TPX = T, P, species
        return self.gas.cp / 1000
    
    def compute(self) -> None:
        molar_heat_capacity: float = np.array([self.get_molar_heat_capacity_fuel({ 'O2': 1, 'N2':3.76 }, T) for T in self.temperature])
        self.enthalpy_fuel: list = FUEL_FORMATION_ENTHALPY + self.temperature * molar_heat_capacity
        molar_heat_capacity: float = np.array([self.get_molar_heat_capacity_fuel({ 'O2': 1, 'N2':3.76 }, T) for T in self.temperature])
        self.enthalpy_oxidier: list = self.temperature * molar_heat_capacity
        print(ct.one_atm)
    def plot(self) -> None:
        graph_1D = Graph_1D()
        graph_1D.setup_axis(
            xlabel='$T-T_0$ [K]', 
            ylabel='Enthalpy [kJ/mol]',
            xmin=TMIN, xmax=TMAX,
            ymin=-1e2, ymax=1e2,
            sci=False
        )
        graph_1D.plot(
            x=self.temperature, 
            y=self.enthalpy_fuel / 1e3,
            label='Fuel',
            color='chartjs_red',
            linestyle='--',
            marker=''
        )
        graph_1D.plot(
            x=self.temperature, 
            y=self.enthalpy_oxidier / 1e3,
            label='Oxidier, Products',
            color='chartjs_purple',
            linestyle='--',
            marker=''
        )
        graph_1D.show(dy=1.1)
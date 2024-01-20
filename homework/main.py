import cantera as ct

from .reaction import Reaction
from .evolution import Evolution
from .enthalpy import Enthalpy

Va0: float = 16.66

def get_fuel() -> dict[float]:
    return { 'C2H6': 1 }

def get_air(e: float) -> dict[float]:
    return { 'O2': (1+e)*Va0/4.76, 'N2': 3.76*(1+e)*Va0/4.76 }

def get_reactants(e: float) -> dict[float]:
    return { **get_fuel(), **get_air(e) }
    
def get_products(e: float) -> dict[float]:
    return { 'CO2': 2, 'H2O': 3, 'N2': 3.76*(1+e)*Va0/4.76, 'O2': e*Va0/4.76 }

def main_1():
    e: float = 0.0
    fuel: dict = get_fuel()
    air: dict = get_air(e)
    reactants: dict = get_reactants(e)
    products: dict = get_products(e)
    
    reaction = Reaction( T0=753, P0=25.12*ct.one_atm, REACTANTS=reactants, PRODUCTS=products )
    
    init_comp: list[float] = reaction.get_composition()
    
    print(f"Molar mass of the fuel: {reaction.get_molar_mass(fuel):.4f} g/mol")
    print(f"Molar mass of the air: {reaction.get_molar_mass(air):.4f} g/mol")
    print(f"Molar mass of the products: {reaction.get_molar_mass(products):.4f} g/mol")
    
    reaction.reactants_reaction()
    
    final_comp: list[float] = reaction.get_composition()
    
    init_comp = [ (k, v) for k, v in sorted(init_comp.items(), key=lambda item: item[1], reverse=True) if v > 1e-3 ]
    final_comp = [ (k, v) for k, v in sorted(final_comp.items(), key=lambda item: item[1], reverse=True) if v > 1e-3 ]
    
    str = ""
    for k, v in init_comp:
        str += f"\n{k:5} {v:.4f}"
    print(f"\nReactants: {str}")
    str = ""
    for k, v in final_comp:
        str += f"\n{k:5} {v:.4f}"
    print(f"\nProducts: {str}")

def main_2():
    evolution = Evolution()
    evolution.compute()
    evolution.plot_concentrations()
    evolution.plot_pressure()
    
def enthalpy():
    enthalpy_sim = Enthalpy()
    enthalpy_sim.compute()
    enthalpy_sim.plot()


def main():
    main_1()
    enthalpy()
    main_2()
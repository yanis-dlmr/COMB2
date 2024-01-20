from dataclasses import dataclass
import cantera as ct

@dataclass
class Reaction:
    T0: float = 298.15
    P0: float = ct.one_atm
    REACTANTS: dict = None
    PRODUCTS: dict = None
    __gas = ct.Solution('gri30.yaml')
    
    def __post_init__(self):
        self.__gas.TPX = self.T0, self.P0, self.REACTANTS

    def get_molar_mass(self, species: dict) -> float:
        """Returns the molar mass of a group of species.
        """
        self.__gas.TPX = self.T0, self.P0, species
        return self.__gas.mean_molecular_weight

    def get_composition(self) -> dict[float]:
        """Returns the composition of the reactants.
        dict = {
            'species name': 'mole fraction'
        }
        """
        return dict(zip(self.__gas.species_names, self.__gas.X))

    def reactants_reaction(self) -> None:
        """Do the reaction of the reactants.
        """
        self.__gas.TPX = self.T0, self.P0, self.REACTANTS
        self.__gas.equilibrate('HP')
        
        print(f"Temperature: {int(self.T0)} -> {int(self.__gas.T)} ({int(self.__gas.T - self.T0)}) K")
        print(f"Pressure: {int(self.P0)} -> {int(self.__gas.P)} ({int(self.__gas.P - self.P0)}) Pa")
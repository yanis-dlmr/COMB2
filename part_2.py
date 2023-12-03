from part_2.assignment import Simulation
from part_2.Gibbs_NO import Gibbs_NO
from part_2.combustion import Combustion

def gibbs():
    sim = Gibbs_NO()
    sim.compute()
    sim.plot()

def combustion():
    sim = Combustion()
    sim.compute_equilibrium_1()
    sim.compute_equilibrium_2()
    sim.compute_evolution()
    sim.plot()

def assigment():
    sim = Simulation()
    sim.compute()
    sim.plot_T_f()
    sim.plot_X_PROD()
    sim.plot_reaction_rate_constant_N1_N2()

if __name__ == '__main__':
    #gibbs()
    #combustion()
    assigment()
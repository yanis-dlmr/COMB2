from part_3.flame_1d_methane import Simulation

COMPUTE: bool = False

def flame_1d_methane() -> None:
    sim = Simulation()
    if COMPUTE:
        sim.compute()
        sim.save_solution()
    else:
        sim.load_solution()
        
    sim.compute_flame_speed('CH4')
    sim.compute_flame_speed('O2')
    sim.compute_flame_speed('N2')
    sim.compute_flame_speed('NO')
    sim.compute_flame_speed('CO2')
    sim.compute_flame_speed('H2O')
    
    sim.plot_mole_fractions()
    sim.plot_flow_velocity()
    sim.plot_flow_density()


if __name__ == '__main__':
    flame_1d_methane()
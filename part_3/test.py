import numpy as np


alphas: np.ndarray = np.linspace(0.0, 1.0, 11)
print(alphas)

for _i, _a in enumerate(alphas):
    
    if (_a == 0 or _a % 0.2 == 0 or _a == 1.0):
        print(_a)
        
        #self.temperature.append(sim.get_temperature())
        #self.z.append(sim.get_grid())
        #self.t_alpha.append(_a)
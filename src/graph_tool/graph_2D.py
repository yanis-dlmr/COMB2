import matplotlib.pyplot as plt
import numpy as np


__all__ = ['Graph_2D']


custom_color = {
    "blue": {
        "normal": "blue",
        "light": "cyan",
    },
    "darkblue": {
        "normal": "darkblue",
        "light": "royalblue",
    },
    "red": {
        "normal": "red",
        "light": "orange",
    },
    "darkred": {
        "normal": "crimson",
        "light": "lightcoral",
    },
    "green": {
        "normal": "green",
        "light": "lightgreen",
    },
    "purple": {
        "normal": "purple",
        "light": "violet",
    },
}


class Graph_2D:
    __tick = 20
    
    def __init__(self):
        self.__fig, __ax = plt.subplots()
        self.__axis = [ __ax ]

    def setup_axis(self, xlabel: str, ylabel: str, xmin: float, xmax: float, ymin: float, ymax: float, axis_number: int = 0, sci: bool = True, tick: int = __tick) -> None:
        """
        Setup the axis of the graph.
        """
        axis = self.__axis[axis_number]
        
        if sci:
            axis.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
            axis.ticklabel_format(style='sci', axis='x', scilimits=(0,0), useMathText=True)
        
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        
        axis.set_xlim(xmin, xmax)
        axis.set_ylim(ymin, ymax)
    
    def plot_map(self, X: list, Y: list, Z: list[list], label: str, axis_number: int = 0, scale_number: int = 7) -> None:
        """
        Plot the graph in 3D.
        """
        plt.style.use('_mpl-gallery-nogrid')
        
        axis = self.__axis[axis_number]
        
        levels = np.linspace(Z.min(), Z.max(), scale_number)

        contourf = axis.contourf(X, Y, Z, levels=levels)

        cbar = self.__fig.colorbar(contourf, ax=axis)
        
    def show(self):
        """
        Show the graph.
        """
        for axis in self.__axis:
            axis.legend()
        plt.grid(True, which="both", linestyle='--', linewidth=0.5, color='grey')
        plt.show()
    
    def add_subplot(self) -> None:
        """
        Add a subplot to the graph.
        """
        __ax = plt.figure().add_subplot()
        self.__axis.append(__ax)
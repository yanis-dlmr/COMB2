import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np


__all__ = ['Graph_3D']


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


class Graph_3D:
    __tick = 20
    
    def __init__(self):
        __ax = plt.figure().add_subplot(projection='3d')
        self.__axis = [ __ax ]

    def setup_axis(self, xlabel: str, ylabel: str, zlabel, xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float, axis_number: int = 0, sci: bool = True, tick: int = __tick) -> None:
        """
        Setup the axis of the graph.
        """
        axis = self.__axis[axis_number]
        
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.set_zlabel(zlabel)
        
        axis.set_xlim(xmin, xmax)
        axis.set_ylim(ymin, ymax)
        axis.set_zlim(zmin, zmax)
        
    def add_subplot(self, axis_number: int = 0) -> None:
        """
        Add a subplot to the graph.
        """
        __ax = plt.figure().add_subplot(projection='3d')
        self.__axis.append(__ax)
        
    
    def plot_3D(self, X: list[list], Y: list[list], Z: list[list], label: str, axis_number: int = 0, rstride: int = 1, cstride: int = 1) -> None:
        """
        Plot the graph in 3D.
        """
        axis = self.__axis[axis_number]

        axis.plot_surface(X, Y, Z, edgecolor='royalblue', lw=0.5, rstride=rstride, cstride=cstride, alpha=0.3)
        
        axis.contour(X, Y, Z, zdir='z', offset=-15, cmap='coolwarm')
        axis.contour(X, Y, Z, zdir='x', offset=-13, cmap='coolwarm')
        axis.contour(X, Y, Z, zdir='y', offset=-13, cmap='coolwarm')
    
    def plot_map(self, X: list[list], Y: list[list], Z: list[list], label: str, axis_number: int = 0) -> None:
        """
        Plot the graph in 3D.
        """
        plt.style.use('_mpl-gallery-nogrid')
        
        axis = self.__axis[axis_number]
        
        X, Y = X[0][:], Y[:,0]
        
        levels = np.linspace(Z.min(), Z.max(), 7)

        axis.contourf(X, Y, Z, levels=levels)
        
    def show(self):
        """
        Show the graph.
        """
        for axis in self.__axis:
            axis.legend()
        plt.grid(True, which="both", linestyle='--', linewidth=0.5, color='grey')
        plt.show()
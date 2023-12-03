import matplotlib.pyplot as plt
import numpy as np


__all__ = ['Graph_1D']


custom_color = {
    "blue": {
        "normal": "blue",
        "light": "cyan",
    },
    "bluesky": {
        "normal": "deepskyblue",
        "light": "lightskyblue",
    },
    "darkblue": {
        "normal": "darkblue",
        "light": "royalblue",
    },
    "red": {
        "normal": "red",
        "light": "orange",
    },
    "lightred": {
        "normal": "lightcoral",
        "light": "lightsalmon",
    },
    "darkred": {
        "normal": "crimson",
        "light": "lightcoral",
    },
    "green": {
        "normal": "green",
        "light": "lightgreen",
    },
    "lightgreen": {
        "normal": "limegreen",
        "light": "lime",
    },
    "purple": {
        "normal": "purple",
        "light": "violet",
    },
    "black": {
        "normal": "black",
        "light": "gray",
    },
    "orange": {
        "normal": "orange",
        "light": "gold",
    },
    "brown": {
        "normal": "brown",
        "light": "sandybrown",
    },
    "pink": {
        "normal": "pink",
        "light": "lightpink",
    },
    "chartjs_blue": {
        "normal": "#36a2eb",
        "light": "#9ad0f5",
    },
    "chartjs_red": {
        "normal": "#ff6384",
        "light": "#ffb1c1",
    },
    "chartjs_green": {
        "normal": "#4bc0c0",
        "light": "#9ad0d0",
    },
    "chartjs_orange": {
        "normal": "#ffcd56",
        "light": "#ffe6b3",
    },
    "chartjs_purple": {
        "normal": "#9966ff",
        "light": "#ccb3ff",
    },
    "chartjs_yellow": {
        "normal": "#ffff00",
        "light": "#ffffb3",
    },
}

class Graph_1D:
    __linestyle = "--"
    __marker = "o"
    __tick = 10
    
    
    def __init__(self):
        self.__fig, __ax = plt.subplots(figsize=(8,6))
        self.__axis = [ __ax ]
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = "Times New Roman"
        plt.rcParams["legend.labelcolor"] = "#363636"
    
    
    def add_axis(self):
        """
        Add an axis to the graph and store it in the graph object.
        """
        ax2 = self.__axis[0].twinx()
        self.__axis.append(ax2)
    
    
    def setup_axis(self, xlabel: str, ylabel: str, xmin: float, xmax: float, ymin: float, ymax: float, axis_number: int = 0, sci: bool = True, tick: int = __tick, color: str = None) -> None:
        """
        Setup the axis of the graph.
        """
        axis = self.__axis[axis_number]
        axis.spines['top'].set_color('#e1e1e1')
        axis.spines['bottom'].set_color('#e1e1e1')
        axis.spines['left'].set_color('#e1e1e1')
        axis.spines['right'].set_color('#e1e1e1')
        axis.set_xlabel(xlabel, fontname="Times New Roman", color="#363636")
        axis.set_ylabel(ylabel, fontname="Times New Roman", color="#363636")
        if sci:
            axis.ticklabel_format(style='sci', axis='x', scilimits=(0,0), useMathText=True)
            axis.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
            for offset in axis.get_yaxis().get_offset_text(), axis.get_xaxis().get_offset_text():
                offset.set_fontname("Times New Roman")
                offset.set_color("#363636")
        else:
            axis.ticklabel_format(style='plain', axis='x', scilimits=(0,0), useMathText=True)
            axis.ticklabel_format(style='plain', axis='y', scilimits=(0,0), useMathText=True)
        axis.set_xlim([xmin, xmax])
        axis.set_ylim([ymin, ymax])
        axis.set_xticks(np.arange(xmin, xmax, (xmax - xmin) / tick))
        axis.set_yticks(np.arange(ymin, ymax, (ymax - ymin) / tick))
        for tick_label in axis.get_xticklabels() + axis.get_yticklabels():
            tick_label.set_fontname("Times New Roman")
            tick_label.set_color("#363636")
        axis.axhline(y=0, color='#e1e1e1', linestyle='-', linewidth=1)
        if color:
            axis.tick_params(axis='y', colors=custom_color[color]["normal"])
    
    def setup_secondary_axis(self, ylabel: str, ymin: float, ymax: float, axis_number: int = 1, sci: bool = True, tick: int = __tick, color: str = None) -> None:
        """
        Setup the axis of the graph.
        """
        axis = self.__axis[axis_number]
        axis.spines['top'].set_color('#e1e1e1')
        axis.spines['bottom'].set_color('#e1e1e1')
        axis.spines['left'].set_color('#e1e1e1')
        axis.spines['right'].set_color('#e1e1e1')
        axis.set_ylabel(ylabel, fontname="Times New Roman", color="#363636")
        if sci:
            axis.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
            offset = axis.get_yaxis().get_offset_text()
            offset.set_fontname("Times New Roman")
            offset.set_color("#363636")
        else:
            axis.ticklabel_format(style='plain', axis='y', scilimits=(0,0), useMathText=True)
        axis.set_ylim([ymin, ymax])
        axis.set_yticks(np.arange(ymin, ymax, (ymax - ymin) / tick))
        for tick_label in axis.get_yticklabels():
            tick_label.set_fontname("Times New Roman")
            tick_label.set_color("#363636")
        axis.axhline(y=0, color='#e1e1e1', linestyle='-', linewidth=1)
        if color:
            axis.tick_params(axis='y', colors=custom_color[color]["normal"])
    
    def plot(self, x: list, y: list, label: str = None, color: str = None, axis_number: int = 0, linestyle: str = __linestyle, marker: str = __marker) -> None:
        """
        Plot the graph.
        """
        axis = self.__axis[axis_number]
        _color = custom_color[color]["normal"] if color else None
        axis.plot(
            x, 
            y, 
            marker=marker, 
            label=label, 
            color=_color,
            linestyle=linestyle
        )
    
    
    def plot_incertitude(self, x: list, y: list, x_incertitude: list, y_incertitude: list, label: str, color: str, axis_number: int = 0) -> None:
        """
        Plot the graph with incertitude.
        """
        if color not in custom_color.keys():
            raise ValueError("The color is not valid.")
        axis = self.__axis[axis_number]
        axis.plot(
            x,
            y,
            linestyle=self.__linestyle,
            marker=self.__marker,
            label=label,
            color=custom_color[color]["normal"]
        )
        axis.errorbar(
            x,
            y,
            xerr=x_incertitude,
            yerr=y_incertitude,
            fmt='none',
            color=custom_color[color]["normal"],
            capsize = 5
        )
        
        
    def plot_interpolation(self, x: list, y: list, interpolation_degree: int, color: str, axis_number: int = 0) -> None:
        """
        Plot the interpolation only, you have to plot the data if you want before or after this function.
        
        Degree 1: linear interpolation y = A * x + B
        Degree 2: quadratic interpolation y = A * x^2 + B * x + C
        """
        if interpolation_degree not in [1, 2]:
            raise ValueError("The interpolation degree must be 1 or 2.")
        
        ax = self.__axis[axis_number]
        
        if interpolation_degree == 1:
            X = np.array(x)
            Y = np.array(y)
            X = X.astype(np.float64)
            Y = Y.astype(np.float64)

            A = np.vstack([X, np.ones(len(X))]).T
            slope, intercept = np.linalg.lstsq(A, Y, rcond=None)[0]
            equation = f'{slope:.5f}x + {intercept:.5f}'
            r_squared = 1 - np.sum((Y - (slope*X + intercept))**2) / ((len(Y) - 1) * np.var(Y, ddof=1))
            
            ax.plot(
                X, 
                slope*X + intercept, 
                color=custom_color[color]["light"],
                linestyle=self.__linestyle, 
                label=f'$Equation: {equation} ; R² = {r_squared:.4f}$'
            )
            
        elif interpolation_degree == 2:
            X = np.array(x)
            Y = np.array(y)
            X = X.astype(np.float64)
            Y = Y.astype(np.float64)
            
            p = np.polyfit(X, Y, 2)
            equation = f'{p[0]:.5f}x^2 + {p[1]:.5f}x + {p[2]:.5f}'
            y_fit = np.polyval(p, X)
            r_squared = 1 - np.sum((Y - y_fit)**2) / ((len(Y) - 1) * np.var(Y, ddof=1))
            
            ax.plot(
                X, 
                y_fit, 
                color=custom_color[color]["light"],
                linestyle=self.__linestyle, 
                label=f'$Equation: {equation} ; R² = {r_squared:.4f}$'
            )
            
            max = -p[1] / (2 * p[0])
            max_value = p[0] * max**2 + p[1] * max + p[2]
            ax.plot(
                max,
                max_value,
                marker='o',
                color=custom_color[color]["light"],
                linestyle='none',
                label=f'$Max: ({max:.4f}, {max_value:.4f})$'
            )
            
        
    def plot_sensitivity(self, x: list, y: list, color: str, label: str, axis_number: int = 0) -> None:
        ax = self.__axis[axis_number]
        
        X = np.array(x)
        Y = np.array(y)
        X = X.astype(np.float64)
        Y = Y.astype(np.float64)
        
        p = np.polyfit(X, Y, 2)
        equation = f'{p[0]:.5f}x^2 + {p[1]:.5f}x + {p[2]:.5f}'
        y_fit = np.polyval(p, X)
        r_squared = 1 - np.sum((Y - y_fit)**2) / ((len(Y) - 1) * np.var(Y, ddof=1))
            
        coeffs = np.polyder(p)
        deriv_equation = f'{coeffs[0]:.5f}x + {coeffs[1]:.5f}'
        y_deriv_fit = np.polyval(coeffs, X)
        ax.plot(
            X, 
            y_deriv_fit, 
            color=custom_color[color]["normal"], 
            linestyle=self.__linestyle, 
            label=f"Equation de la sensibilité pour {label} : ${deriv_equation}$"
        )
        
    def show(self, loc: str = 'upper left', dx: float = 0.3, dy: float = 1.14, ncol: int = 2):
        """
        Show the graph.
        """
        locs = ['upper left', 'upper right']
        for i, axis in enumerate(self.__axis):
            if i == 0:
                j = locs.index(loc)
                axis.legend(bbox_to_anchor=(dx, dy), loc=locs[j], ncol=ncol, frameon=False)
            else:
                j = locs.index(loc)
                j_other = locs.index(loc) - 1
                j_other = (j_other + 1) % 2
                axis.legend(bbox_to_anchor=(dx, dy-0.05), loc=locs[j_other], ncol=ncol, frameon=False)

        for axis in self.__axis:
            axis.grid(True, axis='y', linestyle='-', linewidth=0.5, color='#e1e1e1')
            axis.axhline(y=0, color='#e1e1e1', linestyle='-', linewidth=1)

        self.__axis[0].grid(True, axis='x', linestyle='-', linewidth=0.5, color='#e1e1e1')

        plt.show()

    
    def draw_function(self, function, x_min: float, x_max: float, label: str, color: str, axis_number: int = 0) -> None:
        """
        Draw a function.
        """
        x = np.linspace(x_min, x_max, 1000)
        y = function(x)
        print(x)
        print(y)
        self.plot_curve(x, y, label, color, axis_number)
    
    def plot_curve(self, x: list, y: list, label: str, color: str, axis_number: int = 0, marker: str = "") -> None:
        """
        Plot the graph.
        """
        axis = self.__axis[axis_number]
        axis.plot(
            x, 
            y, 
            label=label,
            color=custom_color[color]["normal"], 
            linestyle=self.__linestyle,
            marker = marker
        )
    
    def add_vertical_line(self, x: float, color: str, label: str, axis_number: int = 0) -> None:
        """
        Add a vertical line.
        """
        axis = self.__axis[axis_number]
        axis.axvline(x=x, color=custom_color[color]["normal"], linestyle=self.__linestyle, label=label)
    
    def plot_area_between_curves(self, curve1: list, curve2: list, color: str, label: str, alpha: float = 0.5, axis_number: int = 0) -> None:
        """
        Plot the area between two curves.
        """
        # if curve1 and curve2 are not None then plot the 2 curves and the area between them
        if curve1 is not None and curve2 is not None:
            axis = self.__axis[axis_number]
            axis.plot(curve1[0], curve1[1], color=custom_color[color]["normal"], linestyle=self.__linestyle)
            axis.plot(curve2[0], curve2[1], color=custom_color[color]["normal"], linestyle=self.__linestyle)
            # fill the area between curve1 and curve2
            axis.fill_between(curve1[0], curve1[1], curve2[1], color=custom_color[color]["light"], alpha=alpha, label=label)
    
    def save(self, filename: str) -> None:
        """
        Save the graph.
        """
        plt.savefig(filename, bbox_inches='tight')
    
    def delete(self) -> None:
        """
        Delete the graph.
        """
        plt.close()
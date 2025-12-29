import numpy as np
import utils
import matplotlib.pyplot as plt

stdRedPotFile : dict[str, float] = {}
BatteryValuesFile : dict[float, float] = {}

# INFO: Parses the .csv files and stores them in dictionaries, that act as 
#       databses. Throws an error if the file is not found.
def parse(name: str) -> dict[str, float]: 
    d: dict[str, float] = {}
    try:
        with open(name, "r") as f:
            next(f)
            for line in f:
                key, value = [s.strip() for s in line.split(",")]
                d[key] = float(value)
            return d
    except OSError:
        print(f'Could not read file: ', name)
        exit()

# INFO: Converts a dictionary with string keys into another dictionary, the 
#       new one using float keys.
def convert(strdict: dict[str, float]) -> dict[float, float]:
    d: dict[float, float] = {}
    for key, value in strdict.items():
        d[float(key)] = value
    return d

# INFO: Parses the number of electrons from the half reaction input, checking 
#       the number before "e-". If no number is found, 1 is returned.
#       new one using float keys.
def getElectrons(elec: str) -> int:
    mult_el = elec[elec.find("e-") - 2]
    if str.isdigit(mult_el):
        n_el = int(mult_el)
    else:
        n_el = 1
    return n_el

# INFO: Parses the stochiometric coefficient from the half reaction input, checking 
#       the first number in the LHS and RHS of the reaction. If no number is found, 
#        1 is returned.
def getStoichCoeffs(elec: str) -> tuple[int, int]:

    parts = elec.split("->")
    ox, red = parts[0].strip(), parts[1].strip()

    if str.isdigit(ox[0]):
        nu_ox = int(ox[0])
    else:
        nu_ox = 1

    if str.isdigit(red[0]):
        nu_red = int(red[0])
    else:
        nu_red = 1

    return nu_red, nu_ox

# INFO: Plots a graph using matplotlibs.pyplot module. Takes as inputs the both 
#       axis of the graph and their labels. The asked input is the name of the 
#       .png that will be saved.
def plotGraph(x: np.ndarray, y: np.ndarray, x_name: str, y_name: str) -> None:
    plt.plot(x, y)
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    name = input("Save the plot as:")
    out = "/tmp/" + name + ".png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot as: {out}')
    plt.close()


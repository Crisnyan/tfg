import math as m
import numpy as np
import scipy.constants as cnt
import utils
import error as e
import matplotlib.pyplot as plt

# INFO: Parses the user input and obtains the standard equilibrium potential via 
#       the standard reduction potential database. Returns both the user input 
#       string and the standard reduction potential of the input half reaction.
def getEo() -> tuple[str, float]:

    try:
        resp = input(f'Input the half reaction happening on the working electrode:\n')
        value = utils.stdRedPotFile.get(resp)
        if value == None:
            resp, value = getEo()
    except Exception as e:
        print("Error:", e)
        resp, value = getEo()
    return resp, value

# INFO: Provides the corresponding equilibrium potential using the nerst equation. 
#       Returns the equilibrium potential and the number of electrons of the reaction.
def getEeq(c_red: float, c_ox: float, T: float) -> tuple[float, int]:
    R = cnt.gas_constant
    F = cnt.value("Faraday constant")

    str_half, Eo = getEo()
    print(f'Eo half-reaction = {Eo}')

    n_el = utils.getElectrons(str_half)
    nu_red, nu_ox = utils.getStoichCoeffs(str_half)

    E_eq =  Eo - (R * T / (n_el * F)) * m.log(c_red ** nu_red / c_ox ** nu_ox)

    return E_eq, n_el


# INFO: Calculates the current density in the range [-abs(overpotential), abs(overpotential)]. 
#       Plots the current density vs the overpotential range obtained.
def ButlerVolmer() -> None:
    F = cnt.value("Faraday constant")
    R = cnt.gas_constant
    T = e.error("temperature")

    c_red = e.error("c_red")
    c_ox = e.error("c_ox")
    j0 = e.error("j0")
    beta_c = e.error("beta_c")
    beta_a = 1.0 - beta_c
    E_ap = e.error("E_ap")
    E_ref = e.error("E_ref")
    E_eq, n_el = getEeq(c_red, c_ox, T)
    print(f'E eq: {E_eq}, n,el: {n_el}')

    eta = E_ap - E_eq - E_ref
    nsteps = 1000 * abs(eta)
    etas = np.linspace(-abs(eta), abs(eta), int(nsteps))
    ja = np.exp(F * n_el * beta_a * etas / (T * R))
    jc = np.exp(F * n_el * -beta_c * etas / (T * R))
    currents = j0 * (ja - jc) 

    utils.plotGraph(etas, currents,'Overpotential (V)', 'Current density (A/m^2)')

import math as m
import numpy as np
import scipy.constants as cnt
import parser
import error as e
import matplotlib.pyplot as plt

def getEeq(c_red: float, c_ox: float, T: float) -> tuple[float, int]:
    R = cnt.gas_constant
    F = cnt.value("Faraday constant")

    str_half, Eo = getEo()
    print(f'Eo half-reaction = {Eo}')

    n_el = parser.getElectrons(str_half)
    nu_red, nu_ox = parser.getStoichCoeffs(str_half)

    E_eq =  Eo - (R * T / (n_el * F)) * m.log(c_red ** nu_red / c_ox ** nu_ox)

    return E_eq, n_el

def getEo() -> tuple[str, float]:

    try:
        resp = input(f'Input the half reaction happening on the working electrode:\n')
        value = parser.stdRedPotFile.get(resp)
        if value == None:
            resp, value = getEo()
    except Exception as e:
        print("Error:", e)
        resp, value = getEo()
    return resp, value

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
    nsteps = 100000 * abs(eta)
    etas = np.linspace(-abs(eta), abs(eta), int(nsteps))
    ja = np.exp(F * n_el * beta_a * etas / (T * R))
    jc = np.exp(F * n_el * -beta_c * etas / (T * R))
    currents = j0 * (ja - jc) 

    plt.plot(etas, currents, 'o')
    plt.ylabel('Current density (A/m^2)')
    plt.xlabel('Overpotential (V)')
    resp = input("Select how the y-axis should be: [lin]ear/[log]arithmic\n")
    name = input("Save the plot as:")
    out = "/tmp/" + name + ".png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot as: {out}')
    plt.close()

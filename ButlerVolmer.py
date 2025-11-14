import math as m
import numpy as np
import scipy.constants as cnt
import parser
import error as e
import matplotlib.pyplot as plt

def getElectrons(elec: str) -> int:
    mult_el = elec[elec.find("e-") - 1]
    print(mult_el)
    if str.isdigit(mult_el):
        n_el = int(mult_el)
    else:
        n_el = 1
    return n_el

def getEeq(c_red_an: float, c_ox_an: float,  c_red_cat: float, c_ox_cat: float, T: float) -> tuple[float, int, int] :
    R = cnt.gas_constant
    F = cnt.value("Faraday constant")
    str_cat, Eo_cat = getVal("E_cat")
    str_an, Eo_an = getVal("E_an")

    n_cat = getElectrons(str_cat)
    n_an = getElectrons(str_an)

    E_cat =  Eo_cat - (R * T / (n_cat * F)) * m.log(c_red_cat / c_ox_cat)
    E_an =  Eo_an - (R * T / (n_an * F)) * m.log(c_red_an / c_ox_an)

    E_cell = E_cat - E_an
    return E_cell, n_cat, n_an 

def getVal(type: str) -> tuple[str, float] :

    if type == "E_cat":
        try:
            resp = input("Input the half reaction of the cathodic electrode:\n")
            value = parser.stdRedPotFile.get(resp)
            if value == None:
                resp, value = getVal("E_cat")
        except Exception as e:
            print("Error:", e)
            resp, value = getVal("E_cat")
    else:
        try:
            resp = input("Input the half reaction of the anodic electrode:\n")
            value = parser.stdRedPotFile.get(resp)
            if value == None:
                resp, value = getVal("E_an")
        except Exception as e:
            print("Error:", e)
            resp, value = getVal("E_an")
    return resp, value

def ButlerVolmer() -> None:
    F = cnt.value("Faraday constant")
    R = cnt.gas_constant
    T = e.error("temperature")

    c_red_an = e.error("c_red_an")
    c_ox_an = e.error("c_ox_an")
    c_red_cat = e.error("c_red_cat")
    c_ox_cat = e.error("c_ox_cat")
    k0 = e.error("k0")
    beta_c = e.error("beta_c")
    beta_a = 1.0 - beta_c
    E_ap = e.error("E_ap")
    E_ref = e.error("E_ref")
    E_eq, n_el_cat, n_el_an = getEeq(c_red_an, c_ox_an, c_red_cat, c_ox_cat, T)
    print(E_eq, n_el_an, n_el_cat)
    nsteps = e.error("n_step")

    eta = E_ap - E_eq - E_ref
    etas = np.linspace(-abs(eta), abs(eta), int(nsteps))
    ja = c_ox_an * np.exp(F * n_el_an * beta_a * etas / (T * R))
    jc = c_red_cat * np.exp(F * n_el_cat * -beta_c * etas / (T * R))
    currents = F * k0 * (ja - jc) 

    plt.plot(etas, currents, 'o')
    plt.ylabel('Current density (A/m^2)')
    plt.xlabel('Overpotential (V)')
    resp = input("Choose how the y-axis should be: [lin]ear/[log]arithmic")
    while resp != "lin" and resp != "log":
        resp = input("The y-axis should be: [lin]ear/[log]arithmic")
    if resp == "log":
        plt.yscale('symlog', linthresh=1e-1)
    name = input("Save the plot as:")
    out = "/tmp/" + name + ".png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot as: {out}')
    plt.close()

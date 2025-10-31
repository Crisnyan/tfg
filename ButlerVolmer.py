import math as m
from numpy.char import isnumeric
import scipy.constants as cnt
import parser
import error as e

def getElectrons(elec: str) -> int:
    mult_el = elec[elec.find("e-") - 1]
    if isnumeric(mult_el):
        n_el = int(mult_el)
    else:
        n_el = 1
    return n_el

def getEeq(c_red: float, c_ox: float, T: float) -> tuple[float, int, int] :
    R = cnt.gas_constant
    F = cnt.value("Faraday constant")
    str_red, Eo_red = getVal("E_red")
    str_ox, Eo_ox = getVal("E_ox")

    n_red = getElectrons(str_red)
    n_ox = getElectrons(str_ox)

    E_red =  Eo_red - (R * T / (n_red * F)) * m.log(c_red / c_ox)
    E_ox =  Eo_ox - (R * T / (n_ox * F)) * m.log(c_red / c_ox)

    E_cell = E_red - E_ox
    return E_cell, n_red, n_ox 

def getVal(type: str) -> tuple[str, float] :

    if type == "E_red":
        try:
            resp = input("Input the half reaction of the cathodic electrode:\n")
            print (resp)
            value = parser.stdRedPotFile[resp]
            if not value:
                resp, value = getVal("E_red")
        except Exception as e:
            print("Error:", e)
            resp, value = getVal("E_red")
    else:
        try:
            resp = input("Input the half reaction of the anodic electrode:\n")
            print (resp)
            value = parser.stdRedPotFile[resp]
            if not value:
                resp, value = getVal("E_ox")
        except Exception as e:
            print("Error:", e)
            resp, value = getVal("E_ox")
    return resp, value

def ButlerVolmer(db: bool) -> float:
    F = cnt.value("Faraday constant")
    R = cnt.gas_constant
    T = e.error("temperature")

    if db:
        c_red_an = e.error("c_red_an")
        c_ox_an = e.error("c_ox_an")
        c_red_cat = e.error("c_red_cat")
        c_ox_cat = e.error("c_ox_cat")
        ko = e.error("ko")
        beta_c = e.error("beta_c")
        beta_a = 1.0 - beta_c
        E_ap_an = e.error("E_ap")
        E_ap_cat = e.error("E_ap")
        E_eq, n_el_ox, n_el_red = getEeq(c_red, c_ox, T)

        eta = E_ap - E_eq
        ja = c_ox * m.exp(F * n_el_ox * beta_a * eta / (T * R))
        jc = c_red * m.exp(F * n_el_red * -beta_c * eta / (T * R))
        return F * ko * (ja - jc)

    else:
        jo = e.error("jo")
        n_el = e.error("n_el")
        eta = e.error("eta")
        beta_c = e.error("beta_c")
        beta_a = 1.0 - beta_c
        ja = m.exp(F * n_el * beta_a * eta / (T * R))
        jc = m.exp(F * n_el * -beta_c * eta / (T * R))
        #print(f"Faraday = {F}, R = {R}, n_el = {n_el}")
        #print(f"eta = {eta}, jo = {jo}")
        #print(f"beta_a = {beta_a}, beta_c = {beta_c}")
        #print(f"ja = {ja}, jc = {jc}")

        return jo * (ja - jc)

import math as m
import scipy.constants as cnt

# INFO: Should get a database of electrodes so there's less inputs needed; the equation works well
def ButlerVolmer() -> float:
    jo = float(input("Set initial j:\n"))
    n_el = int(input("Set the number of electrons involved in the reaction:\n"))
    while n_el < 1:
        print('Error: the number of electrons cannot be less than 1') 
        n_el = int(input("Set the number of electrons involved in the reaction:\n"))
    nu = float(input("Set nu:\n"))
    beta_c = float(input("Set the cathodic symmetry factor:\n"))
    while beta_c < 0 or beta_c > 1:
        print('Error: the symmetry factor must be between 0-1') 
        beta_c = float(input("Set the cathodic symmetry factor:\n"))
    beta_a = 1.0 - beta_c
    temp = float(input("Set the desired temperature:\n"))
    ja = m.exp(cnt.value("Faraday constant") * n_el * beta_a * nu / (temp * cnt.gas_constant))
    jc = m.exp(cnt.value("Faraday constant") * n_el * -beta_c * nu / (temp * cnt.gas_constant))
    # print(f"Faraday = {cnt.value("Faraday constant")}, R = {cnt.gas_constant}, n_el = {n_el}")
    # print(f"nu = {nu}, jo = {jo}")
    # print(f"beta_a = {beta_a}, beta_c = {beta_c}")
    # print(f"ja = {ja}, jc = {jc}")

    return jo * (ja - jc)

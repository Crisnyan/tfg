import numpy as np
import scipy.constants as cnt

def error(type: str) -> int | float:
    temp = 0
    if (type == "n_el"):
        try:
            temp = int(input("Set the number of electrons involved in the reaction:\n"))
            while temp < 1 or temp > 4:
                print('Error: the number of electrons cannot be less than 1 or larger than 4') 
                temp = error("n_el")
        except Exception as e:
            print("Error:", e)
            temp = error("n_el")
    elif (type == "jo"):
        try:
            temp = float(input("Set initial jo:\n"))
        except Exception as e:
            print("Error:", e)
            temp = error("jo")
    elif (type == "nu"):
        try:
            temp = float(input("Set nu:\n"))
        except Exception as e:
            print("Error:", e)
            temp = error("nu")
    elif (type == "beta_c"):
        try:
            temp = float(input("Set the cathodic symmetry factor:\n"))
            while temp < 0 or temp > 1:
                print('Error: the symmetry factor must be between 0-1') 
                temp = error("beta_c")
        except Exception as e:
            print("Error:", e)
            temp = error("nu")
    return temp

def CalculateNuAct() -> float:
    return 0.0
# WARN: Not finished

def CalculateNuConc() -> float:
    return 0.0
# WARN: Not finished

def GetE_Eq() -> float:
# NOTE: extract from CSV all the OCV-SOC
    return 0.0

# WARN: Not finished

def BatteryDischarge() -> float:
    nu_act = CalculateNuAct()
    nu_conc = CalculateNuConc()
    E_eq = GetE_Eq()
    intensity = float(input("Set the intensity:\n"))
    while intensity <= 0:
        print('Error: the number of electrons cannot be less or equal to 0') 
        intensity = float(input("Set the number of electrons involved in the reaction:\n"))
    resistance = float(input("Set the intensity:\n"))
    while resistance <= 0:
        print('Error: the number of electrons cannot be less or equal to 0') 
        intensity = float(input("Set the number of electrons involved in the reaction:\n"))
    return E_eq - nu_act - nu_conc - intensity * resistance

def Battery2():
    F = cnt.value("Faraday constant")
    R = cnt.gas_constant
    T = float(input("Set the desired temperature:\n"))
    n_el = error("n_el")
    E0 = 3.7           # nominal cell potential (V)
    C_nom = 3.0        # capacity in Ah
    R_int = 0.05       # internal resistance (O)
    I = 1.0            # discharge current (A)
    dt = 1.0           # time step (s)
    SOC = 1.0          # start at 100% state-of-charge
    time = 0.0
    voltage = []
    times = []
    while SOC > 0:
        # Nernst-based open-circuit voltage
        E_eq = E0 - (R*T/(n_el*F))*np.log((1.0-SOC)/SOC)
        V = E_eq - I*R_int    # minus ohmic drop (no polarization loss term)
        voltage.append(V)
        times.append(time)
        dSOC = I*dt/(3600.0*C_nom)  # convert to fraction of capacity
        SOC = max(0, SOC - dSOC)
        time += dt
    # Now voltage vs time (or SOC) is in 'voltage' list



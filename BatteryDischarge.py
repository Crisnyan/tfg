import math as m
import scipy.constants as cnt

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

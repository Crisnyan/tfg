import numpy as np
from funcs import ButlerVolmer, BatteryDischarge, LimitIntensity, CyclicVoltamperometry

def ButlerVolmerSelector():
    print("Has seleccionado Butler-Volmer")
    ButlerVolmer()

def BatteryDischargeSelector():
    print("Has seleccionado Descarga de una pila")

def LimitIntensitySelector():
    print("Has seleccionado Intensidad limite")

def CyclicVoltamperometrySelector():
    print("Has seleccionado Voltamperometria ciclica")

def selection():
    print("""Selecciona que tipo de modo quieres usar:
          (1) Butler-Volmer
          (2) Descarga de una pila
          (3) Intensidad limite
          (4) Voltamperometria cilica""")
    resp = int(input())
    return resp

def main():
    resp = selection()
    while resp not in (1, 2, 3, 4):
        print("Valor no permitido, selecciona un numero entre el 1-4")
        resp = selection()
    match resp:
        case 1:
            ButlerVolmerSelector()
        case 2:
            BatteryDischargeSelector()
        case 3:
            LimitIntensitySelector()
        case 4:
            CyclicVoltamperometrySelector()

if __name__ == "__main__":
    main()

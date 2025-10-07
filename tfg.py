from ButlerVolmer import *
from BatteryDischarge import *
from LimitIntensity  import *
from CyclicVoltamperometry import *

def ButlerVolmerSelector():
    print("Butler-Volmer has been selected")
    try:
        print(f"j = {ButlerVolmer()}")
    except Exception as e:
        print("Error:", e)


def BatteryDischargeSelector():
    print("Battery discharge has been selected")
    try:
        print(f"j = {BatteryDischarge()}")
    except Exception as e:
        print("Error:", e)

def LimitIntensitySelector():
    print("Limit intensity has been selected")
    try:
        print(f"j = {LimitIntensity()}")
    except Exception as e:
        print("Error:", e)

def CyclicVoltamperometrySelector():
    print("Cyclic voltamperometry has been selected")
    try:
        print(f"j = {CyclicVoltamperometry()}")
    except Exception as e:
        print("Error:", e)

def selection():
    print("""Selecciona el modo quieras usar:
          (1) Butler-Volmer
          (2) Battery discharge
          (3) Limit intensity
          (4) Cyclic voltamperometry
          (5) Exit""")
    resp = int(input())
    return resp

def main():
    try:
        resp = selection()
        while resp not in (1, 2, 3, 4, 5):
            print("Value not allowed, select a number from 1-4")
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
            case _:
                exit()
    except Exception as e:
        print("Value not allowed, select a number from 1-4")
        main()

if __name__ == "__main__":
    main()

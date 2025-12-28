#!/usr/bin/env python3
from ButlerVolmer import ButlerVolmer
from BatteryDischarge import BatteryDischarge
from RotatingDiskElectrode  import RotatingDiskElectrode
from CyclicVoltammetry import CyclicVoltammetry
import parser

def ButlerVolmerSelector() -> None:
    print("Butler-Volmer has been selected")
    try:
        ButlerVolmer()
    except Exception as e:
        print("Error:", e)

def BatteryDischargeSelector() -> None:
    print("Battery discharge has been selected")
    try:
       BatteryDischarge()
    except Exception as e:
        print("Error:", e)

def RotatingDiskElectrodeSelector() -> None:
    print("Rotating Disk Electrode has been selected")
    try:
        RotatingDiskElectrode()
    except Exception as e:
        print("Error:", e)

def CyclicVoltammetrySelector() -> None:
    print("Cyclic voltamperometry has been selected")
    resp = selection(1)
    while resp not in (1, 2):
        print("Value not allowed, select a valid number")
        resp = selection(1)
    order = 0
    if resp == 2:
        order = selection(2)
        while resp not in (1, 2):
            print("Value not allowed, select a valid number")
            order = selection(2)
    try:
        CyclicVoltammetry(2 * order)
    except Exception as e:
        print("Error:", e)

def alg_selection() -> int:
    resp = int(input())
    return resp

def selection(mode: int) -> int:
    if mode == 0:
        print("""Select the desired use mode:
              (1) Butler-Volmer
              (2) Rotating Disk Electrode
              (3) Battery discharge
              (4) Cyclic voltammetry
              (5) Exit""")
    elif mode == 1:
        print("""Select the desired algorithm:
              (1) Crank-Nicolson + Thomas + Newton
              (2) Runge-Kutta""")
    else:
        print("""Select the desired order:
              (1) Runge-Kutta 2
              (2) Runge-Kutta 4""")
    resp = int(input())
    return resp

def main() -> None:
    try:
        print("Loading", end='')
        parser.stdRedPotFile = parser.parse("standard_potentials.csv")
        print(".", end='')
        parser.BatteryValuesFile = parser.convert(parser.parse("BatteryValues.csv"))
        print(".", end='')
        print(".")
        resp = selection(0)
        while resp not in (1, 2, 3, 4, 5):
            print("Value not allowed, select a number from 1-5")
            resp = selection(0)
        match resp:
            case 1:
                ButlerVolmerSelector()
            case 2:
                RotatingDiskElectrodeSelector()
            case 3:
                BatteryDischargeSelector()
            case 4:
                CyclicVoltammetrySelector()
            case 5:
                exit()
    except Exception as e:
        print("Value not allowed, select a number from 1-5")
        main()

if __name__ == "__main__":
    main()

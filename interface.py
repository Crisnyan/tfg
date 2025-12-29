#!/usr/bin/env python3
from ButlerVolmer import ButlerVolmer
from RotatingDiskElectrode  import RotatingDiskElectrode
from BatteryDischarge import BatteryDischarge
from CyclicVoltammetry import CyclicVoltammetry
import utils

# INFO: Prints the mode selections for the simulations and Runge-Kutta methods.
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

# INFO: Loads the databses and selects which simulation the user desires, handling 
#       parsing errors.
def main() -> None:
    try:
        print("Loading.", end='')
        utils.stdRedPotFile = utils.parse("standard_potentials.csv")
        print(".", end='')
        utils.BatteryValuesFile = utils.convert(utils.parse("BatteryValues.csv"))
        print(".")
        resp = selection(0)
        while resp not in (1, 2, 3, 4, 5):
            print("Value not allowed, select a number from 1-5")
            resp = selection(0)
        match resp:
            case 1:
                print("Butler-Volmer has been selected")
                try:
                    ButlerVolmer()
                except Exception as e:
                    print("Error:", e)
            case 2:
                print("Rotating Disk Electrode has been selected")
                try:
                    RotatingDiskElectrode()
                except Exception as e:
                    print("Error:", e)
            case 3:
                print("Battery discharge has been selected")
                try:
                   BatteryDischarge()
                except Exception as e:
                    print("Error:", e)
            case 4:
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
            case 5:
                exit()
    except Exception as e:
        print("Value not allowed, select a number from 1-5")
        main()

if __name__ == "__main__":
    main()

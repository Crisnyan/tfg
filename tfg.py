from ButlerVolmer import ButlerVolmer
from BatteryDischarge import BatteryDischarge
from LimitIntensity  import LimitIntensity
from CyclicVoltamperometry import CyclicVoltamperometry
import parser

def ButlerVolmerSelector() -> None:
    print("Butler-Volmer has been selected\n")
    resp = input("Use table of Standard reduction potentials?\nyes/no\n")
    while resp != "yes" and resp != "no":
        resp = input("Use table of Standard reduction potentials?\nyes/no\n")
    try:
        if resp == "yes":
            print(f"j = {ButlerVolmer(True)}")
        else:
            print(f"j = {ButlerVolmer(False)}")
    except Exception as e:
        print("Error:", e)


def BatteryDischargeSelector() -> None:
    print("Battery discharge has been selected")
    try:
        print(f"j = {BatteryDischarge()}")
    except Exception as e:
        print("Error:", e)

def LimitIntensitySelector() -> None:
    print("Limit intensity has been selected")
    try:
        print(f"j = {LimitIntensity()}")
    except Exception as e:
        print("Error:", e)

def CyclicVoltamperometrySelector() -> None:
    print("Cyclic voltamperometry has been selected")
    try:
        print(f"j = {CyclicVoltamperometry()}")
    except Exception as e:
        print("Error:", e)

def selection() -> int:
    print("""Select the desired use mode:
          (1) Butler-Volmer
          (2) Battery discharge
          (3) Limit intensity
          (4) Cyclic voltamperometry
          (5) Exit""")
    resp = int(input())
    return resp

def main() -> None:
    try:
        print("Loading", end='')
        parser.stdRedPotFile = parser.parse("standard_potentials.csv")
        print(".", end='')
        parser.BatteryFile = parser.parse("OCVvsSOC.csv")
        parser.BatteryValuesFile = parser.convert(parser.BatteryFile)
        print(".", end='')
        print(".")
        resp = selection()
        while resp not in (1, 2, 3, 4, 5):
            print("Value not allowed, select a number from 1-5")
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
            case 5:
                exit()
    except Exception as e:
        print("Value not allowed, select a number from 1-5")
        main()

if __name__ == "__main__":
    main()

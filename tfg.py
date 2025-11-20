from ButlerVolmer import ButlerVolmer
from BatteryDischarge import BatteryDischarge
from LimitIntensity  import LimitIntensity
from CyclicVoltammetry import CyclicVoltammetry
import parser

def ButlerVolmerSelector() -> None:
    print("Butler-Volmer has been selected\n")
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

def LimitIntensitySelector() -> None:
    print("Limit intensity has been selected")
    try:
        LimitIntensity()
    except Exception as e:
        print("Error:", e)

def CyclicVoltammetrySelector() -> None:
    print("Cyclic voltamperometry has been selected")
    try:
        CyclicVoltammetry()
    except Exception as e:
        print("Error:", e)

def selection() -> int:
    print("""Select the desired use mode:
          (1) Butler-Volmer
          (2) Limit intensity
          (3) Battery discharge
          (4) Cyclic voltammetry
          (5) Exit""")
    resp = int(input())
    return resp

def main() -> None:
    try:
        print("Loading", end='')
        parser.stdRedPotFile = parser.parse("standard_potentials.csv")
        print(".", end='')
        parser.BatteryValuesFile = parser.convert(parser.parse("OCVvsSOC.csv"))
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
                LimitIntensitySelector()
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

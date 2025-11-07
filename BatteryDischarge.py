import numpy as np
import scipy.constants as cnt
import error as e
import matplotlib.pyplot as plt
import parser


def interpOCV(val: float) -> float :
    x1 = 0.0000
    x2 = 0.0000

    print(f'val: {val}')
    keys = list(parser.BatteryValuesFile.keys())
    if val == keys[-1]:
        return parser.BatteryValuesFile[keys[-1]]

    for i in range(1, len(keys)):
        x1 = keys[i - 1]
        x2 = keys[i]
        if x1 <= val <= x2:
            print(f'x1: {x1}, x2: {x2}')
            y1 = parser.BatteryValuesFile[x1]
            y2 = parser.BatteryValuesFile[x2]
            t = (y2 - y1) / (x2 - x1)
            # print(f'x1: {x1}, x2: {x2}, y1: {y1}, y2: {y2}')
            return y1 + t * (val - x1)


def BatteryDischarge() -> None:
    I = e.error("intensity")
    R_int = e.error("resistance")
    Q_nom = e.error("Q_nom")
    State = e.error("SOC") / 100
    t = 0.0
    dt = 1.0
    times = [] 
    voltages = []
    dSOC = - I * dt / (3600.0 * Q_nom)

    print(f'state: {State}, dSOC: {dSOC}')
    while State > 0:
        V = interpOCV(State) - I * R_int
        voltages.append(V)
        times.append(t)
        State = max(0, State + dSOC)
        t += dt
    #print(f'times: {times}')
    #print(f'voltages: {voltages}')
    plt.plot(times, voltages, 'o')
    plt.ylabel('Voltage (V)')
    plt.xlabel('time (s)')
     
    out = "/tmp/battery_discharge.png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot to: {out}')
    plt.close()



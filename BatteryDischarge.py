import numpy as np
import error as e
import matplotlib.pyplot as plt
import parser


def interpOCV(val: float) -> float :
    x1 = 0.0000
    x2 = 0.0000

    keys = list(parser.BatteryValuesFile.keys())
    if val == keys[-1]:
        return parser.BatteryValuesFile[keys[-1]]

    for i in range(1, len(keys)):
        x1 = keys[i - 1]
        x2 = keys[i]
        if x1 <= val <= x2:
            y1 = parser.BatteryValuesFile[x1]
            y2 = parser.BatteryValuesFile[x2]
            t = (y2 - y1) / (x2 - x1)
            return y1 + t * (val - x1)
    return 0.0

def BatteryDischarge() -> None:
    I = e.error("intensity")
    R_int = e.error("resistance")
    Q_nom = e.error("Q_nom")
    State = e.error("SOC") / 100
    endtime = float(e.error("endtime"))
    t = 0.0
    dt = 1.0
    voltages = []
    dSOC = - I * dt / (Q_nom)
    if endtime == 0.0:
        endtime = 9e34

    # print(f'state: {State}, dSOC: {dSOC}')
    while State > 0 and t <= endtime:
        V = interpOCV(State) - I * R_int
        voltages.append(V)
        State += dSOC
        t += dt
    voltages = np.array(voltages)
    times = np.arange(len(voltages) * dt)
    #print(f'times: {times}')
    #print(f'voltages: {voltages}')
    plt.plot(times, voltages, 'o')
    plt.ylabel('Voltage (V)')
    plt.xlabel('time (s)')
    name = input("Save the plot as:")
    out = "/tmp/" + name + ".png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot as: {out}')
    plt.close()

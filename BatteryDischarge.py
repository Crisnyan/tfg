import numpy as np
import scipy.constants as cnt
import error as e
import matplotlib.pyplot as plt
import parser


def interpOCV(val: float) -> float :
    step = 10
    val *= 100
    floor = val // step * step
    ceil = floor + step
    pond = (step - val % step) / step

    n1 = parser.BatteryFile[str(floor)]
    n2 = parser.BatteryFile[str(ceil)]

    return n1 * pond + n2 * (1 - pond)

def BatteryDischarge() -> float:
    I = e.error("intensity")
    R_int = e.error("resistance")
    Q_nom = e.error("Q_nom")
    State = e.error("SOC") / 100
    t = 0.0
    dt = 1.0
    times = [] 
    voltages = []
    dSOC = - I * dt / (3600.0 * Q_nom)

    while State > 0:
        V = interpOCV(State) - I * R_int
        voltages.append(V)
        times.append(t)
        State = max(0, State + dSOC)
        t += dt
    plt.plot(times, voltages, 'o')
    plt.ylabel('Voltage (V)')
    plt.xlabel('time (s)')

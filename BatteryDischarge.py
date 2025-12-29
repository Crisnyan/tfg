import numpy as np
import error as e
import matplotlib.pyplot as plt
import utils


# INFO: Linearly interpolates the input SOC's corresponding potential value, 
#       using as interpolation points experimentally  obtained values for a 
#       Li-ion battery from the BatteryValues database.  Returns the corresponding 
#       linearly interpolated potential.
def interpOCV(val: float) -> float :
    x1 = 0.0000
    x2 = 0.0000

    keys = list(utils.BatteryValuesFile.keys())
    if val == keys[-1]:
        return utils.BatteryValuesFile[keys[-1]]

    for i in range(1, len(keys)):
        x1 = keys[i - 1]
        x2 = keys[i]
        if x1 <= val <= x2:
            y1 = utils.BatteryValuesFile[x1]
            y2 = utils.BatteryValuesFile[x2]
            t = (y2 - y1) / (x2 - x1)
            return y1 + t * (val - x1)
    return 0.0

# INFO: Calculates the terminal voltage of the Li-ion battery, looping until
#       the time from endtime is reached or the battery's fully discharged.
#       Plots the discharge curve (Voltage vs time):
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

    while State > 0 and t <= endtime:
        V = interpOCV(State) - I * R_int
        voltages.append(V)
        State += dSOC
        t += dt
    voltages = np.array(voltages)
    times = np.arange(len(voltages) * dt)

    utils.plotGraph(times, voltages, 'time (s)', 'Voltage (V)')

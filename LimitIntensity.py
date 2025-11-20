import numpy as np
import error as e
import scipy.constants as cnt
import matplotlib.pyplot as plt

def LimitIntensity():
    n_el = e.error("n_el")
    F = cnt.value("Faraday constant")
    D_coef = e.error("diff")
    C0 = e.error("C0")
    nu = e.error("viscosity")
    A = e.error("area")
    rpm_max = e.error("rpm_max")

    rpm = np.linspace(100, rpm_max, 50)
    omega = 2 * np.pi * rpm / 60.0

    #Levich 
    delta = 1.61 * (D_coef**(1/3)) * (omega**(-0.5)) * (nu**(1/6))
    i_lim = n_el * F * D_coef * C0 * A / delta

    plt.plot(rpm, i_lim)
    plt.xlabel("RPM")
    plt.ylabel("Limit current (A)")
    name = input("Save the plot as:")
    out = "/tmp/" + name + ".png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot as: {out}')
    plt.close()
    print(f"With {rpm[-1]:.0f} RPM, the limit current is {i_lim[-1]:.4f} A")

import scipy.constants as cnt
import error as e

def LimitIntensity() -> float:
    # given parameters (SI)
    n_el = e.error("n_el")
    F = cnt.value("Faraday constant")
    D = 1e-9          # m^2/s
    c0 = 1.0          # mol/m^3
    delta = e.error("delta")

    # Planar steady-state limiting current (if you assume a fixed delta)
    return n_el * F * D * c0 / delta

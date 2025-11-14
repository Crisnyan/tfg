import numpy as np
import math as m
import scipy.constants as cnt
import error as e

def init_time_potential(E_start: float, E_vertex: float, v: float, dt: float, n_cycles: int) -> tuple[np.ndarray, np.ndarray]:
    E_min = E_start - E_vertex
    E_max = E_start + E_vertex

    t_cycle = 4.0 * E_vertex / v 
    t_max = n_cycles * t_cycle
    t = np.linspace(0.0, t_max, int(round(t_max / dt)) + 1)

    phase_shift = (E_start - E_min) / v
    u = np.mod(t + phase_shift, t_cycle) / t_cycle
    tri = 1.0 - 2.0 * np.abs(u - 0.5)
    E = E_min + tri * (E_max - E_min)

    return t, E

def CyclicVoltammetry():
    F = cnt.value("Faraday constant")
    R = cnt.gas_constant
    T = e.error("temperature")
    nx = 300

    E_start = e.error("E_start")
    E_vertex = abs(e.error("E_vertex") - E_start)
    srate = e.error("srate")
    n_cycles = e.error("n_cycles")
    bulk_Co = e.error("Cob")
    bulk_Cr = e.error("Crb")
    D_coef = e.error("diff")
    k0 = e.error("k0")
    beta_c = e.error("beta_c")
    beta_a = 1.0 - beta_c
    n_el = e.error("n_el")

    t_max = 4.0 * n_cycles * E_vertex / srate
    x_max = 3.0 * m.sqrt(2.0 * D_coef * t_max)
    dx = x_max / nx
    dt = dx**2 / (5.0 * D_coef)
    t, E = init_time_potential(E_start, E_vertex, srate, dt, int(n_cycles))

    Co = np.full(nx + 1, bulk_Co)
    Cr = np.full(nx + 1, bulk_Cr)

    j = np.empty_like(t)
    #print(f't_max: {t_max}\n, x_max: {x_max}\n, dt: {dt}\n, dx: {dx}\n\n')
    #print(f'time array: {t}\n, time len: {len(t)}\n,  potential array: {E}\n,  potential len: {len(E)}\n, Co array: {Co}\n, Cr array: {Cr}\n')
    for i in range(len(t) - 1):
        E_next = E[i + 1]  # applied potential at next time (we update concentrations to t+dt under E_next)
        eta = E_next - E0  # overpotential (V) = applied - formal potential

        # surface concentrations at current time (before update)
        cO_s = Co[0]
        cR_s = Cr[0]

        # --- Butler-Volmer-like surface reaction rate (mol m^-2 s^-1) ---
        # using didactic form: rate = k0 * ( c_O * exp(-alpha*nF*eta/RT) - c_R * exp((1-alpha)*nF*eta/RT) )
        # k0 units [m/s], concentrations [mol/m^3] => rate [mol/m^2/s]
        exp_f = m.exp(-beta_a * n_el * F * eta / (R * T))   # forward factor (reduction)
        exp_b = m.exp(beta_c * n_el * F * eta / (R * T))    # backward factor (oxidation)
        # note: signs chosen so that positive r corresponds to net reduction O -> R
        r_surf = k0 * (cO_s * exp_f - cR_s * exp_b)  # mol / (m^2 s)

        # convert to current density (A / m^2)
        j[i + 1] = n_el * F * r_surf

        # --- explicit FD update for the concentrations ---
        # copy arrays to build next-time-step values
        Co_new = Co.copy()
        Cr_new = Cr.copy()

        # interior nodes update (i = 1 .. N-2) using FTCS
        # c_i^{n+1} = c_i^n + coef * (c_{i+1}^n - 2 c_i^n + c_{i-1}^n)
        Co_new[1:-1] = Co[1:-1] + coef * (Co[2:] - 2.0 * Co[1:-1] + Co[:-2])
        Cr_new[1:-1] = Cr[1:-1] + coef * (Cr[2:] - 2.0 * Cr[1:-1] + Cr[:-2])

        # far-field boundary condition (rightmost node) : Dirichlet to bulk concentrations
        Co_new[-1] = bulk_Co
        Cr_new[-1] = bulk_Cr

        # surface node update (i = 0) using ghost-point / flux condition given by r_surf:
        # For species O:
        # Flux of O at surface J_O = - r_surf  (mol / m^2 s) because r positive consumes O -> R
        # one-sided second derivative approximation (see derivation in comments earlier):
        # d2c/dx2|_0 ~ 2*(c1 - c0)/dx^2 + 2*J_O/(D*dx)
        J_O = - r_surf   # mol / (m^2 s)
        d2cO_dx2_0 = 2.0 * (Co[1] - Co[0]) / (dx * dx) + 2.0 * J_O / (D_coef * dx)
        Co_new[0] = Co[0] + D_coef * dt * d2cO_dx2_0

        # For species R:
        # Flux of R at surface J_R = + r_surf
        J_R = + r_surf
        d2cR_dx2_0 = 2.0 * (Cr[1] - Cr[0]) / (dx * dx) + 2.0 * J_R / (D_coef * dx)
        Cr_new[0] = Cr[0] + D_coef * dt * d2cR_dx2_0

        # numerical clamping to avoid tiny negative values due to roundoff
        Co_new = np.maximum(Co_new, 0.0)
        Cr_new = np.maximum(Cr_new, 0.0)

        # commit update for next iteration
        Co = Co_new
        Cr = Cr_new


    #    j[i] = n_el * F * r


    





    #r = k0 * ( Cob * m.exp(-beta_c * n_el * F * eta / (R * T)) - Crb * m.exp((beta_a) * n_el * F * eta / (R * T)) )



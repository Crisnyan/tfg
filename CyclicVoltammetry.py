import numpy as np
import math as m
import scipy.constants as cnt
import error as e
import parser
import matplotlib.pyplot as plt


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

def getVal() -> tuple[int, float] :

    value = None
    n_el = 0
    while value == None:
        resp = input("Input the half reaction corresponding to the oxidation:\n")
        value = parser.stdRedPotFile.get(resp)
        n_el = parser.getElectrons(resp)
    return n_el, value

def thomas(a_temp, b_temp, c_temp, d_temp):
    n = len(d_temp)                                      # number of unknowns

    a, b, c, d = map(np.array, (a_temp, b_temp, c_temp, d_temp)) 
    for i in range(1, n):                           # FORWARD ELIMINATION
        m = a[i] / b[i - 1]                           # compute multiplier that eliminates a[i]
        b[i] = b[i] - m * c[i - 1]                    # update pivot in row i
        d[i] = d[i] - m * d[i - 1]                    # update RHS in row i

    x = np.zeros(n)
    x[-1] = d[-1] / b[-1]                           # back-sub: solve last equation

    for i in range(n - 2, -1, -1):                    # BACK SUBSTITUTION
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]        # solve each earlier x[i]
    return x

def rhs(conc: np.ndarray, c0: float, bulk: float, lam: float) -> np.ndarray:
    res = 0.5 * lam * conc[0:-2] + (1.0 - lam) * conc[1:-1] + 0.5 * lam * conc[2:]
    # incorporate new surface value Co0 into first RHS entry (ghost/elimination)
    res[0] += 0.5 * lam * c0
    # incorporate far-field Dirichlet (node nx) — use bulk values at new time (assumed equal to bulk)
    res[-1] += 0.5 * lam * bulk
    return res

def reaction_rate_and_partials(Co0: float, Cr0: float, E_ap: float, params: dict[str, int | float]):

    R = params["R"]
    T = params["T"]
    F = params["F"]
    k0 = params["k0"]
    n_el = params["n_el"]
    beta_c = params["beta_c"]
    beta_a = 1 - beta_c
    E_std = params["E_std"]

    f = F / (R * T)
    xi = f * (E_ap - E_std)

    k_red = k0 * m.exp(-beta_c * n_el * xi)
    k_ox  = k0 * m.exp((beta_a) * n_el * xi)

    r = k_ox * Cr0 - k_red * Co0

    # analytic partials (derived in the notes)
    # dr/dCo0 = (1 - alpha) * k0 * (j_a + j_c) / Co0
    # dr/dCr0 = alpha * r / Cr0
    dr_dCo0 = -k_red
    dr_dCr0 = k_ox
    # protect division by zero

    return r, dr_dCo0, dr_dCr0


def newton_thomas(Co: np.ndarray, Cr: np.ndarray, E_ap: float, lam: float, params: dict[str, int | float]):
    tol=1e-9
    max_iter=12

    dx = params["dx"]
    D_coef = params["D"]

    # initial guess: previous surface concentrations
    n_unknowns = Co.size - 2

    low_diag = np.full(n_unknowns, -0.5 * lam) 
    mid_diag = np.full(n_unknowns, (1.0 + lam))
    upp_diag = np.full(n_unknowns, -0.5 * lam) 

    u = np.array([Co[0], Cr[0]])
    for it in range(max_iter):
        Co0, Cr0 = u[0], u[1]

        do = rhs(Co, Co0, Co[-1], lam)
        dr = rhs(Cr, Cr0, Cr[-1], lam)
        # solve linear systems for interiors
        Co_interior = thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), do)
        Cr_interior = thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), dr)

        # compute residuals R_O and R_R (ghost-point derived)
        # R_O = Co0 - Co_prev0 - lam*( 2*(C1 - Co0) - (2*dx/D)*r_surf )

        r_surf, dr_dCo0, dr_dCr0 = reaction_rate_and_partials(Co0, Cr0, E_ap, params)

        Co1 = Co_interior[0]
        Cr1 = Cr_interior[0]

        #resid_ox = Co0 - Co[0] - lam * ( 2.0 * (Co1 - Co0) - (2.0 * dx / D_coef) * r_surf )
        #resid_red = Cr0 - Cr[0] - lam * ( 2.0 * (Cr1 - Cr0) + (2.0 * dx / D_coef) * r_surf )
        resid_ox = Co0 - Co1 - dx / D_coef * r_surf
        resid_red = Cr0 - Cr1 + dx / D_coef * r_surf # Signo opuesto para reducción


        if max(abs(resid_ox), abs(resid_red)) < tol:
            # assemble full arrays for return: insert interior into full profile
            return Co0, Cr0, Co_interior, Cr_interior, r_surf

        b_sens = np.zeros(n_unknowns)
        b_sens[0] = 0.5 * lam
        s_vec = thomas(low_diag, mid_diag, upp_diag, b_sens)
        s_factor = s_vec[0] # dC1 / dC0
        
        # J = [ dR_ox/dCo0   dR_ox/dCr0 ]
        #     [ dR_red/dCo0  dR_red/dCr0 ]
        
        # R_ox = Co0 - Co1(Co0) - scale * r(Co0, Cr0)
        # dR_ox/dCo0 = 1 - s_factor - scale * dr/dCo0
        # dR_ox/dCr0 = - scale * dr/dCr0
        
        J11 = 1.0 - s_factor - dx / D_coef * dr_dCo0
        J12 = - dx / D_coef * dr_dCr0
        
        # R_red = Cr0 - Cr1(Cr0) + scale * r(Co0, Cr0)
        # dR_red/dCo0 = scale * dr/dCo0
        # dR_red/dCr0 = 1 - s_factor + scale * dr/dCr0
        
        J21 = dx / D_coef * dr_dCo0
        J22 = 1.0 - s_factor + dx / D_coef * dr_dCr0
        
        J = np.array([[J11, J12], [J21, J22]])
        R_vec = np.array([resid_ox, resid_red])
        
        try:
            delta = np.linalg.solve(J, -R_vec)
        except np.linalg.LinAlgError:
            delta = -R_vec * 0.1 # Fallback
            
        # --- CORRECCIÓN CRÍTICA: ACTUALIZAR U ---
        u += delta
        
        # Evitar valores negativos físicos
        u[u < 0] = 1e-15

    return u[0], u[1], Co_interior, Cr_interior, r_surf
        # # 5) compute sensitivity scalars s_C = dC1/dCo0 and s_R = dR1/dCr0
        # #    Solve A * s = b_C where b_C has 0.5*lam at index0 and zeros elsewhere.
        # bC = np.zeros(n_unknowns, dtype=float)
        # bC[0] = 0.5 * lam
        # sC_vec = np.array(thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), bC.copy()), dtype=float)
        # s_C = sC_vec[0]
        #
        # bR = np.zeros(n_unknowns, dtype=float)
        # bR[0] = 0.5 * lam
        # sR_vec = np.array(thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), bR.copy()), dtype=float)
        # s_R = sR_vec[0]
        #
        # # 6) build analytic Jacobian (2x2) using derived formulas
        # pref = lam * (2.0 * dx / D_coef)
        # J11 = 1.0 - lam * (2.0 * s_C - 2.0) + pref * dr_dCo0
        # J12 = pref * dr_dCr0
        # J21 = - pref * dr_dCo0
        # J22 = 1.0 - lam * (2.0 * s_R - 2.0) - pref * dr_dCr0
        #
        # J = np.array([[J11, J12], [J21, J22]], dtype=float)
        # Rvec = np.array([resid_ox, resid_red], dtype=float)
        #
        # # 7) solve small 2x2 linear system J du = -R
        # #    fallback to small step if solve fails
        # try:
        #     du = np.linalg.solve(J, -Rvec)
        # except np.linalg.LinAlgError:
        #     # singular Jacobian: take small damped step towards negative residual
        #     du = -Rvec * 1e-3
    # return Co0, Cr0, Co_interior, Cr_interior # after max_iter: return best effort
    #Co0, Cr0 = max(u[0], eps), max(u[1], eps)

    # recompute final interiors for returned surface values
    # d_o = 0.5 * lam * Co[0:-2] + (1.0 - lam) * Co[1:-1] + 0.5 * lam * Co[2:]
    # d_r = 0.5 * lam * Cr[0:-2] + (1.0 - lam) * Cr[1:-1] + 0.5 * lam * Cr[2:]
    # d_o[0] += 0.5 * lam * Co0
    # d_r[0] += 0.5 * lam * Cr0
    # d_o[-1] += 0.5 * lam * Co[-1]
    # d_r[-1] += 0.5 * lam * Cr[-1]
    #
    # Co_int = np.array(thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), d_o.copy()), dtype=float)
    # Cr_int = np.array(thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), d_r.copy()), dtype=float)
    #
    # r_final, _, _ = reaction_rate_and_partials(Co0, Cr0, E_ap, params)
    #return Co0, Cr0, Co_int, Cr_int, r_final

def CyclicVoltammetry():
    F = cnt.value("Faraday constant")
    R = cnt.gas_constant
    T = 298.15
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
    n_el, E_std = getVal()

    t_max = 4.0 * n_cycles * E_vertex / srate
    x_max = 3.0 * m.sqrt(2.0 * D_coef * t_max)
    dx = x_max / nx
    dt = dx**2 / (5.0 * D_coef)
    t, E = init_time_potential(E_start, E_vertex, srate, dt, int(n_cycles))

    Co = np.full(nx + 1, bulk_Co)
    Cr = np.full(nx + 1, bulk_Cr)
    lam =  D_coef * dt / (dx**2)

    params = {
        "k0": k0,
        "beta_c": beta_c,      # your beta_a
        "n_el": n_el,
        "F": F,
        "R": R,
        "T": T,
        "dx": dx,
        "D": D_coef,
        "E_std": E_std
    }

    j = np.empty_like(t)
    #print(f't_max: {t_max}\n, x_max: {x_max}\n, dt: {dt}\n, dx: {dx}\n\n')
    #print(f'time array: {t}\n, time len: {len(t)}\n,  potential array: {E}\n,  potential len: {len(E)}\n, Co array: {Co}\n, Cr array: {Cr}\n')
    for i in range(len(t) - 1):

        newCo = np.empty_like(Co)
        newCr = np.empty_like(Cr)

        Co0_new, Cr0_new, Co_int, Cr_int, r_surf = newton_thomas(Co, Cr, E[i], lam, params)
        
        # Reconstruir vectores completos para el siguiente paso
        newCo = np.zeros_like(Co)
        newCr = np.zeros_like(Cr)
        
        newCo[0] = Co0_new
        newCo[1:-1] = Co_int
        newCo[-1] = bulk_Co # Dirichlet
        
        newCr[0] = Cr0_new
        newCr[1:-1] = Cr_int
        newCr[-1] = bulk_Cr # Dirichlet

        Co = newCo
        Cr = newCr

        j[i] = n_el * F * r_surf


    j[-1] = j[-2]
    plt.plot(E, j)
    plt.xlabel('Potential (V)')
    plt.ylabel('Current density (A/m^2)')
    name = input("Save the plot as:")
    out = "/tmp/" + name + ".png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot as: {out}')
    plt.close()

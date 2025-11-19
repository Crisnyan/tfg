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

def thomas(a, b, c, d):
    n = len(d)                                      # number of unknowns

    for i in range(1, n):                           # FORWARD ELIMINATION
        w = a[i] / b[i - 1]                           # compute multiplier that eliminates a[i]
        b[i] = b[i] - w * c[i - 1]                    # update pivot in row i
        d[i] = d[i] - w * d[i - 1]                    # update RHS in row i
    x = [0]*n                                       # allocate solution vector
    x[-1] = d[-1] / b[-1]                           # back-sub: solve last equation
    for i in range(n - 2, -1, -1):                    # BACK SUBSTITUTION
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]        # solve each earlier x[i]
    return x

# def rhs(Co: np.ndarray, Cr: np.ndarray, u: np.ndarray, Cob: float, Crb: float, lam: float) -> tuple[np.ndarray, np.ndarray]:
#     do = 0.5 * lam * Co[0:2] + (1.0 - lam) * Co[1:-1] + 0.5 * lam * Co[2:]
#     dr = 0.5 * lam * Cr[0:2] + (1.0 - lam) * Cr[1:-1] + 0.5 * lam * Cr[2:]
#     # incorporate new surface value Co0 into first RHS entry (ghost/elimination)
#     do[0] += 0.5 * lam * u[0]
#     dr[0] += 0.5 * lam * u[1]
#     # incorporate far-field Dirichlet (node nx) — use bulk values at new time (assumed equal to bulk)
#     do[-1] += 0.5 * lam * Cob
#     dr[-1] += 0.5 * lam * Crb
#     return do, dr

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
    n_el = params["n_el"]
    beta_c = params["beta_c"]
    beta_a = 1 - beta_c
    E_std = params["E_std"]
    k0 = params["k0"]

    eps = 1e-20
    Co0 = max(Co0, eps)
    Cr0 = max(Cr0, eps)

    B = (n_el * F) / (R * T)
    Delta = B * (E_ap - E_std)

    # Gamma = B*eta = B*(E_ap - E_std) + ln(Cr0/Co0)
    Gamma = Delta + m.log(Cr0 / Co0)

    # microscopic forward/backward "factors"
    j_a = Co0 * m.exp(beta_a * Gamma)
    j_c = Cr0 * m.exp(-(1.0 - beta_a) * Gamma)   # note (1-beta_a) == beta_c

    r = k0 * (j_a - j_c)

    # analytic partials (derived in the notes)
    # dr/dCo0 = (1 - alpha) * k0 * (j_a + j_c) / Co0
    # dr/dCr0 = alpha * r / Cr0
    dr_dCo0 = (1.0 - beta_a) * k0 * (j_a + j_c) / Co0
    # protect division by zero
    dr_dCr0 = beta_a * (r / Cr0)

    return r, dr_dCo0, dr_dCr0


def newton_thomas(Co: np.ndarray, Cr: np.ndarray, E_ap: float, lam: float, params: dict[str, int | float]):
    tol=1e-9
    max_iter=12
    eps = 1e-20

    R = params["R"]
    T = params["T"]
    E_std = params["E_std"]
    dx = params["dx"]
    D_coef = params["D"]
    F = params["F"]
    n_el = params["n_el"]
    beta_c = params["beta_c"]

    # initial guess: previous surface concentrations
    u = np.array([Co[0], Cr[0]])
    n_unknowns = Co.size - 2

    low_diag = np.full(n_unknowns, -0.5 * lam) 
    mid_diag = np.full(n_unknowns, (1 + lam))
    upp_diag = np.full(n_unknowns, -0.5 * lam) 

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

        resid_ox = Co0 - Co[0] - lam * ( 2.0 * (Co1 - Co0) - (2.0 * dx / D_coef) * r_surf )
        resid_red = Cr0 - Cr[0] - lam * ( 2.0 * (Cr1 - Cr0) + (2.0 * dx / D_coef) * r_surf )


        if max(abs(resid_ox), abs(resid_red)) < tol:
            # assemble full arrays for return: insert interior into full profile
            return Co0, Cr0, Co_interior, Cr_interior, r_surf

        # 5) compute sensitivity scalars s_C = dC1/dCo0 and s_R = dR1/dCr0
        #    Solve A * s = b_C where b_C has 0.5*lam at index0 and zeros elsewhere.
        bC = np.zeros(n_unknowns, dtype=float)
        bC[0] = 0.5 * lam
        sC_vec = np.array(thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), bC.copy()), dtype=float)
        s_C = sC_vec[0]

        bR = np.zeros(n_unknowns, dtype=float)
        bR[0] = 0.5 * lam
        sR_vec = np.array(thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), bR.copy()), dtype=float)
        s_R = sR_vec[0]

        # 6) build analytic Jacobian (2x2) using derived formulas
        pref = lam * (2.0 * dx / D_coef)
        J11 = 1.0 - lam * (2.0 * s_C - 2.0) + pref * dr_dCo0
        J12 = pref * dr_dCr0
        J21 = - pref * dr_dCo0
        J22 = 1.0 - lam * (2.0 * s_R - 2.0) - pref * dr_dCr0

        J = np.array([[J11, J12], [J21, J22]], dtype=float)
        Rvec = np.array([resid_ox, resid_red], dtype=float)

        # 7) solve small 2x2 linear system J du = -R
        #    fallback to small step if solve fails
        try:
            du = np.linalg.solve(J, -Rvec)
        except np.linalg.LinAlgError:
            # singular Jacobian: take small damped step towards negative residual
            du = -Rvec * 1e-3

        # # finite-difference Jacobian (2x2) - perturb each component
        # J = np.zeros((2,2), dtype=float)
        # R_base = np.array([R_O, R_R], dtype=float)
        # for col in range(2):
        #     u_pert = u.copy()
        #     perturb = max(1e-8, 1e-6 * abs(u[col]) )
        #     u_pert[col] += perturb
        #     Co0_p, Cr0_p = u_pert[0], u_pert[1]
        #     r_surf_p = reaction_rate(Co0_p, Cr0_p, E_ap)
        #     # build RHS with perturbed surface
        #     do = rhs(Co, Co0, Co[-1], lam)
        #     dr = rhs(Cr, Cr0, Cr[-1], lam)
        #
        #     Co_interior = thomas(np.full(do.size, -0.5 * lam) , np.full(do.size, (1 + lam)), np.full(do.size, -0.5 * lam), do)
        #     Cr_interior = thomas(np.full(dr.size, -0.5 * lam) , np.full(dr.size, (1 + lam)), np.full(dr.size, -0.5 * lam), dr)
        #
        #     C1p = Co_int_p[0]
        #     CR1p = Cr_int_p[0]
        #     R_O_p = Co0_p - Co_prev_arr[0] - lam * ( 2.0*(C1p - Co0_p) - (2.0 * dx / D_coef) * r_surf_p )
        #     R_R_p = Cr0_p - Cr_prev_arr[0] - lam * ( 2.0*(CR1p - Cr0_p) + (2.0 * dx / D_coef) * r_surf_p )
        #     R_pert = np.array([R_O_p, R_R_p], dtype=float)
        #     J[:, col] = (R_pert - R_base) / perturb
        #
        # # solve small linear system J * du = -R to update u
        # try:
        #     du = np.linalg.solve(J, -R_base)
        # except np.linalg.LinAlgError:
        #     # Jacobian singular — very small step or fallback
        #     du = -R_base * 1e-3
        #
        # # optional damping: ensure residual decreases
        # lambda_d = 1.0
        # accepted = False
        # for _ in range(5):
        #     u_trial = u + lambda_d * du
        #     Co0_t, Cr0_t = u_trial[0], u_trial[1]
        #     r_surf_t = reaction_rate(Co0_t, Cr0_t, E_ap)
        #     # quick residual eval (cheap)
        #     dO_t = np.empty(n_unknown, dtype=float)
        #     dR_t = np.empty(n_unknown, dtype=float)
        #     for j in range(n_unknown):
        #         i_node = j + 1
        #         dO_t[j] = 0.5*lam * Co_prev_arr[i_node - 1] + (1.0 - lam) * Co_prev_arr[i_node] + 0.5*lam * Co_prev_arr[i_node + 1]
        #         dR_t[j] = 0.5*lam * Cr_prev_arr[i_node - 1] + (1.0 - lam) * Cr_prev_arr[i_node] + 0.5*lam * Cr_prev_arr[i_node + 1]
        #     dO_t[0] += 0.5 * lam * Co0_t
        #     dR_t[0] += 0.5 * lam * Cr0_t
        #     dO_t[-1] += 0.5 * lam * bulk_Co
        #     dR_t[-1] += 0.5 * lam * bulk_Cr
        #     Co_int_t = thomas(a_template, b_template, c_template, dO_t)
        #     Cr_int_t = thomas(a_template, b_template, c_template, dR_t)
        #     C1t = Co_int_t[0]
        #     CR1t = Cr_int_t[0]
        #     R_O_t = Co0_t - Co_prev_arr[0] - lam * ( 2.0*(C1t - Co0_t) - (2.0 * dx / D_coef) * r_surf_t )
        #     R_R_t = Cr0_t - Cr_prev_arr[0] - lam * ( 2.0*(CR1t - Cr0_t) + (2.0 * dx / D_coef) * r_surf_t )
        #     if max(abs(R_O_t), abs(R_R_t)) < res_norm:
        #         accepted = True
        #         u = u_trial
        #         break
        #     lambda_d *= 0.5
        #
        # if not accepted:
        #     # If damping failed to reduce residual, do a small step
        #     u += 0.1 * du
        #
    # end Newton iterations
    # # if we get here without returning, accept last iterate (best effort)
    # Co0, Cr0 = u[0], u[1]
    # # recompute interiors for last iterate
    # dO_last = np.empty(n_unknown, dtype=float)
    # dR_last = np.empty(n_unknown, dtype=float)
    # for j in range(n_unknown):
    #     i_node = j + 1
    #     d_O_last[j] = 0.5*lam * Co_prev[i_node - 1] + (1.0 - lam) * Co_prev[i_node] + 0.5*lam * Co_prev[i_node + 1]
    #     d_R_last[j] = 0.5*lam * Cr_prev[i_node - 1] + (1.0 - lam) * Cr_prev[i_node] + 0.5*lam * Cr_prev[i_node + 1]
    # d_O_last[0] += 0.5 * lam * Co0
    # d_R_last[0] += 0.5 * lam * Cr0
    # d_O_last[-1] += 0.5 * lam * bulk_Co
    # d_R_last[-1] += 0.5 * lam * bulk_Cr
    # Co_interior = thomas(a_template, b_template, c_template, d_O_last)
    # Cr_interior = thomas(a_template, b_template, c_template, d_R_last)
    # return Co0, Cr0, Co_interior, Cr_interior # after max_iter: return best effort
    Co0, Cr0 = max(u[0], eps), max(u[1], eps)

    # recompute final interiors for returned surface values
    d_o = 0.5 * lam * Co[0:-2] + (1.0 - lam) * Co[1:-1] + 0.5 * lam * Co[2:]
    d_r = 0.5 * lam * Cr[0:-2] + (1.0 - lam) * Cr[1:-1] + 0.5 * lam * Cr[2:]
    d_o[0] += 0.5 * lam * Co0
    d_r[0] += 0.5 * lam * Cr0
    d_o[-1] += 0.5 * lam * Co[-1]
    d_r[-1] += 0.5 * lam * Cr[-1]

    Co_int = np.array(thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), d_o.copy()), dtype=float)
    Cr_int = np.array(thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), d_r.copy()), dtype=float)

    r_final, _, _ = reaction_rate_and_partials(Co0, Cr0, E_ap, params)
    return Co0, Cr0, Co_int, Cr_int, r_final

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
    n_el, E_std = getVal()
    print(f'n_el {n_el}, E_std {E_std}')

    t_max = 4.0 * n_cycles * E_vertex / srate
    x_max = 3.0 * m.sqrt(2.0 * D_coef * t_max)
    dx = x_max / nx
    dt = dx**2 / (5.0 * D_coef)
    t, E = init_time_potential(E_start, E_vertex, srate, dt, int(n_cycles))
    print("pasa")

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

    print("setea parametros")
    j = np.empty_like(t)
    #print(f't_max: {t_max}\n, x_max: {x_max}\n, dt: {dt}\n, dx: {dx}\n\n')
    #print(f'time array: {t}\n, time len: {len(t)}\n,  potential array: {E}\n,  potential len: {len(E)}\n, Co array: {Co}\n, Cr array: {Cr}\n')
    for i in range(len(t) - 1):

        newCo = np.empty_like(Co)
        newCr = np.empty_like(Cr)

        newCo[0], newCr[0], newCo[1:-1], newCr[1:-1], r_surf = newton_thomas(Co, Cr, E[i], lam, params)

        newCo[-1] = bulk_Co
        newCr[-1] = bulk_Cr

        Co = newCo
        Cr = newCr

        j[i] = n_el * F * r_surf


    plt.plot(E, j)
    plt.xlabel('Potential (V)')
    plt.ylabel('Current density (A/m^2)')
    name = input("Save the plot as:")
    out = "/tmp/" + name + ".png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot as: {out}')
    plt.close()

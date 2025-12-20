import numpy as np
import math as m
import scipy.constants as cnt
import error as e
import parser
import matplotlib.pyplot as plt
import matplotlib.animation as anim

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
    n = len(d_temp)

    a, b, c, d = map(np.array, (a_temp, b_temp, c_temp, d_temp)) 
    for i in range(1, n):
        m = a[i] / b[i - 1]
        b[i] = b[i] - m * c[i - 1]
        d[i] = d[i] - m * d[i - 1]

    x = np.zeros(n)
    x[-1] = d[-1] / b[-1]

    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    return x

def rhs(conc: np.ndarray, c0: float, bulk: float, lam: float) -> np.ndarray:
    res = 0.5 * lam * conc[0:-2] + (1.0 - lam) * conc[1:-1] + 0.5 * lam * conc[2:]
    res[0] += 0.5 * lam * c0
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

    dr_dCo0 = -k_red
    dr_dCr0 = k_ox

    return r, dr_dCo0, dr_dCr0

def newton_thomas(Co: np.ndarray, Cr: np.ndarray, E_ap: float, lam: float, params: dict[str, int | float]):
    tol=1e-9
    max_iter=12

    dx = params["dx"]
    D_coef = params["D"]

    n_unknowns = Co.size - 2

    low_diag = np.full(n_unknowns, -0.5 * lam) 
    mid_diag = np.full(n_unknowns, (1.0 + lam))
    upp_diag = np.full(n_unknowns, -0.5 * lam) 

    u = np.array([Co[0], Cr[0]])
    for it in range(max_iter):
        Co0, Cr0 = u[0], u[1]

        do = rhs(Co, Co0, Co[-1], lam)
        dr = rhs(Cr, Cr0, Cr[-1], lam)

        Co_interior = thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), do)
        Cr_interior = thomas(low_diag.copy(), mid_diag.copy(), upp_diag.copy(), dr)

        r_surf, dr_dCo0, dr_dCr0 = reaction_rate_and_partials(Co0, Cr0, E_ap, params)

        Co1 = Co_interior[0]
        Cr1 = Cr_interior[0]

        resid_ox = Co0 - Co1 - dx / D_coef * r_surf
        resid_red = Cr0 - Cr1 + dx / D_coef * r_surf # INFO: Signo opuesto para reducción


        if max(abs(resid_ox), abs(resid_red)) < tol:
            return Co0, Cr0, Co_interior, Cr_interior, r_surf

        b_sens = np.zeros(n_unknowns)
        b_sens[0] = 0.5 * lam
        s_vec = thomas(low_diag, mid_diag, upp_diag, b_sens)
        s_factor = s_vec[0]

        J11 = 1.0 - s_factor - dx / D_coef * dr_dCo0
        J12 = - dx / D_coef * dr_dCr0

        J21 = dx / D_coef * dr_dCo0
        J22 = 1.0 - s_factor + dx / D_coef * dr_dCr0

        J = np.array([[J11, J12], [J21, J22]])
        R_vec = np.array([resid_ox, resid_red])

        try:
            delta = np.linalg.solve(J, -R_vec)
        except np.linalg.LinAlgError:
            delta = -R_vec * 0.1

        u += delta
        u[u < 0] = 1e-15

    return u[0], u[1], Co_interior, Cr_interior, r_surf

def solve_surface_analytic(Co1: float, Cr1: float, E_ap: float, params: dict[str, int | float]):
    """
    Resuelve analíticamente las concentraciones en superficie (nodo 0)
    basándose en el flujo desde el nodo 1 y las cinéticas Butler-Volmer.
    Sistema: D(C1 - C0)/dx = Rate
    """
    R = params["R"]
    T = params["T"]
    F = params["F"]
    k0 = params["k0"]
    n_el =  params["n_el"]
    beta_c = params["beta_c"]
    beta_a = 1 - beta_c
    E_std = params["E_std"]
    D_coef = params["D"]
    dx = params["dx"]

    f = F / (R * T)
    xi = f * (E_ap - E_std)
    k_red = k0 * m.exp(-beta_c * n_el * xi)
    k_ox  = k0 * m.exp(beta_a * n_el * xi)

    # INFO: Coeficiente difusivo 'alpha' = D/dx
    alpha = D_coef / dx

    # INFO: Sistema lineal de 2 ecuaciones:
    # INFO: 1) alpha*(Co1 - Co0) = k_ox*Cr0 - k_red*Co0  => (alpha + k_red)*Co0 - k_ox*Cr0 = alpha*Co1
    # INFO: 2) alpha*(Cr1 - Cr0) = -(k_ox*Cr0 - k_red*Co0) => -k_red*Co0 + (alpha + k_ox)*Cr0 = alpha*Cr1
    
    A = np.array([
        [alpha + k_red, -k_ox],
        [-k_red, alpha + k_ox]
    ])
    b_vec = np.array([alpha * Co1, alpha * Cr1])
    
    # INFO: Resolver para [Co0, Cr0]
    res = np.linalg.solve(A, b_vec)
    return res[0], res[1], k_ox, k_red

def get_derivatives(t: float, Co_int: np.ndarray, Cr_int: np.ndarray, E_now: float, params: dict[str, int | float]) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Calcula dC/dt para los nodos internos usando el método de líneas.
    Co_int, Cr_int: Arrays de nodos internos (excluyendo 0 y N)
    E_func_params: tuple (E_start, E_vertex, v, etc) para calcular E al tiempo t
    """
    D = params["D"]
    dx = params["dx"]
    bulk_Co = params["bulk_Co"]
    bulk_Cr = params["bulk_Cr"]

    Co0, Cr0, k_ox, k_red = solve_surface_analytic(Co_int[0], Cr_int[0], E_now, params)

    Co_full = np.concatenate(([Co0], Co_int, [bulk_Co]))
    Cr_full = np.concatenate(([Cr0], Cr_int, [bulk_Cr]))

    # INFO: 4. Calcular Laplaciano (Difusión) dC/dt = D * d2C/dx2
    # INFO: Vectorizado: (C[i+1] - 2C[i] + C[i-1]) / dx^2
    dCo_dt = D * (Co_full[2:] - 2 * Co_full[1:-1] + Co_full[:-2]) / (dx ** 2)
    dCr_dt = D * (Cr_full[2:] - 2 * Cr_full[1:-1] + Cr_full[:-2]) / (dx ** 2)

    # INFO: Calculamos la corriente actual (flux) para guardarla externamente si fuera necesario,
    # INFO: pero aquí solo devolvemos derivadas de estado.
    r_surf = k_ox * Cr0 - k_red * Co0

    return dCo_dt, dCr_dt, r_surf

def CyclicVoltammetry(order: int):
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

    t_max = 4.0 * n_cycles * E_vertex / srate
    x_max = 3.0 * m.sqrt(2.0 * D_coef * t_max)
    dx = x_max / nx
    dt = dx ** 2 / (5.0 * D_coef)
    t, E = init_time_potential(E_start, E_vertex, srate, dt, int(n_cycles))
    x = np.arange(0, nx + 1, 1)

    Co = np.full(nx + 1, bulk_Co)
    Cr = np.full(nx + 1, bulk_Cr)
    j = np.empty_like(t)

    Co_anim = np.full((len(t), nx + 1), bulk_Co)
    Cr_anim = np.full((len(t), nx + 1), bulk_Cr)

    lam =  D_coef * dt / (dx ** 2)

    params = {
        "k0": k0,
        "beta_c": beta_c,
        "n_el": n_el,
        "F": F,
        "R": R,
        "T": T,
        "dx": dx,
        "D": D_coef,
        "E_std": E_std,
        "bulk_Co": bulk_Co,
        "bulk_Cr": bulk_Cr
    }

    dt2 = dt / 2
    dt6 = dt / 6

    if (order == 0):
        for i in range(len(t) - 1):

            Co0_new, Cr0_new, Co_int, Cr_int, r_surf = newton_thomas(Co, Cr, E[i], lam, params)


            newCo = np.concatenate(([Co0_new], Co_int, [bulk_Co]))
            newCr = np.concatenate(([Cr0_new], Cr_int, [bulk_Cr]))

            Co = newCo
            Cr = newCr

            j[i] = n_el * F * r_surf
            Co_anim[i + 1] = newCo
            Cr_anim[i + 1] = newCr

    elif order == 2:

        Co_int = Co[1:-1]
        Cr_int = Cr[1:-1]

        for i in range(len(t) - 1):

            k1_Co, k1_Cr, r_surf = get_derivatives(t[i], Co_int, Cr_int, E[i], params)

            Co_pred = Co_int + dt * k1_Co
            Cr_pred = Cr_int + dt * k1_Cr

            k2_Co, k2_Cr, _ = get_derivatives(t[i + 1], Co_pred, Cr_pred, E[i + 1], params)

            Co_int = Co_int + dt2 * (k1_Co + k2_Co)
            Cr_int = Cr_int + dt2 * (k1_Cr + k2_Cr)

            Co0_new, Cr0_new, _, _ = solve_surface_analytic(Co_int[0], Cr_int[0], E[i + 1], params)

            Co_anim[i + 1] = np.concatenate(([Co0_new], Co_int, [bulk_Co]))
            Cr_anim[i + 1] = np.concatenate(([Cr0_new], Cr_int, [bulk_Cr]))

            j[i] = n_el * F * r_surf

    elif order == 4:

        Co_int = Co[1:-1]
        Cr_int = Cr[1:-1]

        for i in range(len(t) - 1):

            k1_Co, k1_Cr, r_surf = get_derivatives(t[i], Co_int, Cr_int, E[i], params)
            k2_Co, k2_Cr, _ = get_derivatives(t[i] + dt2, Co_int + dt2 * k1_Co, Cr_int + dt2 * k1_Cr, 0.5 * (E[i] + E[i + 1]), params)
            k3_Co, k3_Cr, _ = get_derivatives(t[i] + dt2, Co_int + dt2 * k2_Co, Cr_int + dt2 * k2_Cr, 0.5 * (E[i] + E[i + 1]), params)
            k4_Co, k4_Cr, _ = get_derivatives(t[i + 1], Co_int + dt * k3_Co, Cr_int + dt * k3_Cr, E[i + 1], params)

            Co_int = Co_int + dt6 * (k1_Co + 2 * k2_Co + 2 * k3_Co + k4_Co)
            Cr_int = Cr_int + dt6 * (k1_Cr + 2 * k2_Cr + 2 * k3_Cr + k4_Cr)

            Co0_new, Cr0_new, _, _ = solve_surface_analytic(Co_int[0], Cr_int[0], E[i + 1], params)

            Co_anim[i + 1] = np.concatenate(([Co0_new], Co_int, [bulk_Co]))
            Cr_anim[i + 1] = np.concatenate(([Cr0_new], Cr_int, [bulk_Cr]))

            j[i] = n_el * F * r_surf

        j[-1] = j[-2]
    parser.save_plot(E, j, 'Potential (V)', 'Current density (A/m^2)')
    #animate_diffusion(x, Co_anim, Cr_anim)

def animate_diffusion(x, Co_anim, Cr_anim):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    
    max_conc = max(np.max(Co_anim), np.max(Cr_anim)) * 1.1
    
    ax1.set_xlabel("distance (m)")
    ax1.set_ylabel("Conc. Oxidized (mol/m^3)")
    ax1.set_xlim(0, x[-1])     
    ax1.set_ylim(0, max_conc)   
    
    line1, = ax1.plot([], [], 'b-', lw=2) 
    
    ax2.set_xlabel("distance (m)")
    ax2.set_ylabel("Conc. Reduced (mol/m^3)")
    ax2.set_xlim(0, x[-1])
    ax2.set_ylim(0, max_conc)
    
    line2, = ax2.plot([], [], 'r-', lw=2)

    time_text = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)

    def update(frame):
        
        line1.set_data(x, Co_anim[frame])
        
        line2.set_data(x, Cr_anim[frame])
        
        time_text.set_text(f'Time Step: {frame}')
        
        return line1, line2, time_text

    ani = anim.FuncAnimation(
        fig, 
        update, 
        frames=len(Co_anim), 
        interval=20, 
        blit=True
    )

    plt.tight_layout()
    
    plt.show()



    #plt.figure()
    #plt.subplot(1, 2, 1)
    #plt.plot(x, Co_anim)
    #plt.xlabel("distance (m)")
    #plt.ylabel("concentration of the oxidized species (mol/m^3)")

    #plt.subplot(1, 2, 2)
    #plt.plot(x, Cr_anim)
    #plt.xlabel("distance (m)")
    #plt.ylabel("concentration of the reduced species (mol/m^3)")


    #name = input("Save the plot as:")
    #out = "/tmp/" + name + ".png"
    #plt.savefig(out, dpi=150, bbox_inches='tight')
    #print(f'Saved plot as: {out}')
    #plt.close()
    #parser.save_plot(E, Co_surf, 'Potential (V)', 'Concentration of the oxidized species (mol/m^3)')
    #parser.save_plot(E, Cr_surf, 'Potential (V)', 'Concentration of the reduced species (mol/m^3)')

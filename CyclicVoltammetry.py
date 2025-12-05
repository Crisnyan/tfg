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

def getVal() -> tuple[int, float]:

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
        m_factor = a[i] / b[i - 1]
        b[i] = b[i] - m_factor * c[i - 1]
        d[i] = d[i] - m_factor * d[i - 1]

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
        resid_red = Cr0 - Cr1 + dx / D_coef * r_surf 
        
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

def build_new(first: float, middle: np.ndarray, last: float) -> np.ndarray:
    new = np.empty(len(middle) + 2)

    new[0] = first
    new[1:-1] = middle
    new[-1] = last

    return new

# --- NUEVAS FUNCIONES PARA RK2 / RK4 ---

def solve_surface_analytic(Co1, Cr1, E_ap, params):
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

    # Coeficiente difusivo 'alpha' = D/dx
    alpha = D_coef / dx

    # Sistema lineal de 2 ecuaciones:
    # 1) alpha*(Co1 - Co0) = k_ox*Cr0 - k_red*Co0  => (alpha + k_red)*Co0 - k_ox*Cr0 = alpha*Co1
    # 2) alpha*(Cr1 - Cr0) = -(k_ox*Cr0 - k_red*Co0) => -k_red*Co0 + (alpha + k_ox)*Cr0 = alpha*Cr1
    
    A = np.array([
        [alpha + k_red, -k_ox],
        [-k_red, alpha + k_ox]
    ])
    b_vec = np.array([alpha * Co1, alpha * Cr1])
    
    # Resolver para [Co0, Cr0]
    res = np.linalg.solve(A, b_vec)
    return res[0], res[1], k_ox, k_red

def get_derivatives(t, Co_int, Cr_int, E_func_params, params):
    """
    Calcula dC/dt para los nodos internos usando el método de líneas.
    Co_int, Cr_int: Arrays de nodos internos (excluyendo 0 y N)
    E_func_params: tuple (E_start, E_vertex, v, etc) para calcular E al tiempo t
    """
    # 1. Recuperar E al tiempo t
    E_start, E_vertex, v, phase_shift, t_cycle, E_min, E_max = E_func_params
    
    # Calcular E actual (lógica duplicada de init_time_potential para un t arbitrario)
    u = np.mod(t + phase_shift, t_cycle) / t_cycle
    tri = 1.0 - 2.0 * np.abs(u - 0.5)
    E_now = E_min + tri * (E_max - E_min)

    # 2. Recuperar Constantes
    D = params["D"]
    dx = params["dx"]
    
    # 3. Resolver Frontera Superficial (BC)
    # Necesitamos Co[0] y Cr[0]. Usamos Co[1] que es Co_int[0]
    Co0, Cr0, k_ox, k_red = solve_surface_analytic(Co_int[0], Cr_int[0], E_now, params)
    
    # 4. Construir vectores completos temporales para calcular laplaciano
    # BC bulk es constante
    bulk_Co = params["bulk_Co"]
    bulk_Cr = params["bulk_Cr"]
    
    # Concatenar: [Superficie, Internos..., Bulk]
    Co_full = np.concatenate(([Co0], Co_int, [bulk_Co]))
    Cr_full = np.concatenate(([Cr0], Cr_int, [bulk_Cr]))
    
    # 5. Calcular Laplaciano (Difusión) dC/dt = D * d2C/dx2
    # Vectorizado: (C[i+1] - 2C[i] + C[i-1]) / dx^2
    dCo_dt = D * (Co_full[2:] - 2*Co_full[1:-1] + Co_full[:-2]) / (dx**2)
    dCr_dt = D * (Cr_full[2:] - 2*Cr_full[1:-1] + Cr_full[:-2]) / (dx**2)
    
    # Calculamos la corriente actual (flux) para guardarla externamente si fuera necesario,
    # pero aquí solo devolvemos derivadas de estado.
    r_surf = k_ox * Cr0 - k_red * Co0
    
    return dCo_dt, dCr_dt, r_surf

def CyclicVoltammetry(order: int):
    # --- SETUP INICIAL (Igual que antes) ---
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
    dt = dx ** 2 / (5.0 * D_coef) # Condición de estabilidad para explícitos
    
    t, E = init_time_potential(E_start, E_vertex, srate, dt, int(n_cycles))

    # Parámetros para pasar a funciones
    params = {
        "k0": k0, "beta_c": beta_c, "n_el": n_el, "F": F, "R": R, "T": T,
        "dx": dx, "D": D_coef, "E_std": E_std, 
        "bulk_Co": bulk_Co, "bulk_Cr": bulk_Cr
    }

    # Arrays de estado completo
    Co = np.full(nx + 1, bulk_Co)
    Cr = np.full(nx + 1, bulk_Cr)
    j = np.empty_like(t)

    # --- LÓGICA DE SOLVER ---
    
    if order == 0:
        # --- CRANK-NICOLSON (IMPLÍCITO) ---
        print("Using Crank-Nicolson (Implicit)...")
        lam = D_coef * dt / (dx**2)
        for i in range(len(t) - 1):
            Co0_new, Cr0_new, Co_int, Cr_int, r_surf = newton_thomas(Co, Cr, E[i+1], lam, params) # E[i+1] porque es implícito backward
            
            Co = build_new(Co0_new, Co_int, bulk_Co)
            Cr = build_new(Cr0_new, Cr_int, bulk_Cr)
            j[i] = n_el * F * r_surf
            
        j[-1] = j[-2]

    else:
        # --- MÉTODOS EXPLÍCITOS (RK) ---
        print(f"Using Runge-Kutta Order {order}...")
        
        # Preparamos parámetros para reconstruir E dentro de RK
        E_min = E_start - E_vertex
        E_max = E_start + E_vertex
        t_cycle = 4.0 * E_vertex / srate
        phase_shift = (E_start - E_min) / srate
        E_params = (E_start, E_vertex, srate, phase_shift, t_cycle, E_min, E_max)

        # Trabajamos solo con los nodos internos para la integración ODE
        # Los nodos de borde se recalculan dentro de get_derivatives
        Co_int = Co[1:-1]
        Cr_int = Cr[1:-1]

        for i in range(len(t) - 1):
            t_curr = t[i]
            
            if order == 2:
                # RK2 (Heun / Punto medio)
                k1_Co, k1_Cr, r_surf = get_derivatives(t_curr, Co_int, Cr_int, E_params, params)
                
                # Estado predictor
                Co_pred = Co_int + dt * k1_Co
                Cr_pred = Cr_int + dt * k1_Cr
                
                k2_Co, k2_Cr, _ = get_derivatives(t_curr + dt, Co_pred, Cr_pred, E_params, params)
                
                Co_int = Co_int + (dt / 2.0) * (k1_Co + k2_Co)
                Cr_int = Cr_int + (dt / 2.0) * (k1_Cr + k2_Cr)
                
                # Guardamos corriente (basada en el estado al inicio del paso o promedio)
                j[i] = n_el * F * r_surf

            elif order == 4:
                # RK4
                # k1
                k1_Co, k1_Cr, r_surf = get_derivatives(t_curr, Co_int, Cr_int, E_params, params)
                
                # k2
                k2_Co, k2_Cr, _ = get_derivatives(t_curr + 0.5*dt, Co_int + 0.5*dt*k1_Co, Cr_int + 0.5*dt*k1_Cr, E_params, params)
                
                # k3
                k3_Co, k3_Cr, _ = get_derivatives(t_curr + 0.5*dt, Co_int + 0.5*dt*k2_Co, Cr_int + 0.5*dt*k2_Cr, E_params, params)
                
                # k4
                k4_Co, k4_Cr, _ = get_derivatives(t_curr + dt, Co_int + dt*k3_Co, Cr_int + dt*k3_Cr, E_params, params)
                
                # Update
                Co_int = Co_int + (dt / 6.0) * (k1_Co + 2*k2_Co + 2*k3_Co + k4_Co)
                Cr_int = Cr_int + (dt / 6.0) * (k1_Cr + 2*k2_Cr + 2*k3_Cr + k4_Cr)
                
                j[i] = n_el * F * r_surf

        j[-1] = j[-2]

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(E, j)
    # Invertir eje x (convención IUPAC a veces, o US convención) - opcional según preferencia
    plt.xlabel('Potential (V)')
    plt.ylabel('Current density (A/m^2)')
    plt.title(f'Cyclic Voltammetry (Order={order})')
    plt.grid(True, alpha=0.3)
    
    name = input("Save the plot as: ")
    out = "/tmp/" + name + ".png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved plot as: {out}')
    plt.close()

# Para ejecutar, llama a la función con el orden deseado, ej:
# CyclicVoltammetry(0) # Crank-Nicolson
# CyclicVoltammetry(2) # RK2
# CyclicVoltammetry(4) # RK4

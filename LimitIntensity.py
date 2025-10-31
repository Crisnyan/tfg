import scipy.constants as cnt

def error(type: str) -> int | float:
    temp = 0
    if (type == "n_el"):
        try:
            temp = int(input("Set the number of electrons involved in the reaction:\n"))
            while temp < 1 or temp > 4:
                print('Error: the number of electrons cannot be less than 1 or larger than 4') 
                temp = error("n_el")
        except Exception as e:
            print("Error:", e)
            temp = error("n_el")
    elif (type == "jo"):
        try:
            temp = float(input("Set initial jo:\n"))
        except Exception as e:
            print("Error:", e)
            temp = error("jo")
    elif (type == "area"):
        try:
            temp = float(input("Set area of the electrode (m^2):\n"))
            while temp < 0:
                print('Error: Area cannot be negative') 
            while temp > 1:
                print('Error: Area must be less than 1 m^2') 
        except Exception as e:
            print("Error:", e)
            temp = error("area")
    elif (type == "delta"):
        try:
            temp = float(input("Set diffusion layer:\n"))
            while temp < 0:
                print('Error: The diffusion layer cannot be negative') 
            while temp > 1:
                print('Error: The diffusion layer must be less than 1') 
        except Exception as e:
            print("Error:", e)
            temp = error("delta")
    elif (type == "beta_c"):
        try:
            temp = float(input("Set the cathodic symmetry factor:\n"))
            while temp < 0 or temp > 1:
                print('Error: the symmetry factor must be between 0-1') 
                temp = error("beta_c")
        except Exception as e:
            print("Error:", e)
            temp = error("nu")
    return temp

def LimitIntensity() -> float:
    # given parameters (SI)
    n_el = error("n_el")
    F = cnt.value("Faraday constant")
    area = error("area")
    D = 1e-9          # m^2/s
    c0 = 1.0          # mol/m^3
    delta = error("delta")

    # Planar steady-state limiting current (if you assume a fixed delta)
    return n_el * F * area * D * c0 / delta

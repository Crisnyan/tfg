def error(type: str) -> int | float:
    temp = 0
    match type:
        case "n_el":
            try:
                temp = int(input("Set the number of electrons involved in the reaction:\n"))
                while temp < 1 or temp > 4:
                    print('Error: the number of electrons cannot be less than 1 or larger than 4') 
                    temp = error("n_el")
            except Exception as e:
                print("Error:", e)
                temp = error("n_el")
        case "temperature":
            try:
                temp = float(input("Set the desired temperature:\n"))
                if temp == 0:
                    temp = 1e-12
                while temp < 0:
                    print('Error: the temperature cannot be negative') 
                    temp = error("temperature")
            except Exception as e:
                print("Error:", e)
                temp = error("temperature")
        case "c_red_an":
            try:
                temp = float(input("Set the concentration of the reduced species in the anode (1 for solids):\n"))
                if temp == 0:
                    temp = 1e-12
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("c_red_an")
            except Exception as e:
                print("Error:", e)
                temp = error("c_red_an")
        case "c_ox_an":
            try:
                temp = float(input("Set the concentration of the oxidized species in the anode (1 for solids):\n"))
                if temp == 0:
                    temp = 1e-12
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("c_ox_an")
            except Exception as e:
                print("Error:", e)
                temp = error("c_ox_an")
        case "c_red_cat":
            try:
                temp = float(input("Set the concentration of the reduced species in the cathode (1 for solids):\n"))
                if temp == 0:
                    temp = 1e-12
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("c_red_cat")
            except Exception as e:
                print("Error:", e)
                temp = error("c_red_cat")
        case "c_ox_cat":
            try:
                temp = float(input("Set the concentration of the oxidized species in the cathode (1 for solids):\n"))
                if temp == 0:
                    temp = 1e-12
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("c_ox_cat")
            except Exception as e:
                print("Error:", e)
                temp = error("c_ox_cat")
        case "intensity":
            try:
                temp = float(input("Set the intensity of the battery:\n"))
                if temp == 0:
                    temp = 1e-12
                if temp < 0:
                    temp *= -1
            except Exception as e:
                print("Error:", e)
                temp = error("intensity")
        case "resistance":
            try:
                temp = float(input("Set the internal resistance of the battery:\n"))
                if temp == 0:
                    temp = 1e-12
                while temp < 0:
                    print('Error: resistance cannot be negative') 
                    temp = error("resistance")
            except Exception as e:
                print("Error:", e)
                temp = error("resistance")
        case "Q_nom":
            try:
                temp = float(input("Set the nominal capacity of the battery (mAh):\n"))
                if temp == 0:
                    temp = 1e-12
                while temp < 0:
                    print('Error: the nominal capacity cannot be negative') 
                    temp = error("Q_nom")
            except Exception as e:
                print("Error:", e)
                temp = error("Q_nom")
        case "jo":
            try:
                temp = float(input("Set initial jo:\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("jo")
        case "beta_c":
            try:
                temp = float(input("Set the cathodic symmetry factor:\n"))
                while temp < 0 or temp > 1:
                    print('Error: the symmetry factor must be between 0-1') 
                    temp = error("beta_c")
            except Exception as e:
                print("Error:", e)
                temp = error("beta_c")
        case "delta":
            try:
                temp = float(input("Set diffusion layer:\n"))
                while temp < 0:
                    print('Error: The diffusion layer cannot be negative') 
                while temp > 1:
                    print('Error: The diffusion layer must be 1 or less') 
            except Exception as e:
                print("Error:", e)
                temp = error("delta")
        case "SOC":
            try:
                temp = float(input("Set state of charge of the battery (%):\n"))
                while temp < 0:
                    print('Error: The state of charge of the battery cannot be negative') 
                while temp > 100:
                    print('Error: The state of charge must be 100% or less') 
            except Exception as e:
                print("Error:", e)
                temp = error("SOC")
        case "ko":
            try:
                temp = float(input("Set ko (mol/m^3):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("ko")
        case "eta":
            try:
                temp = float(input("Set the overpotential (V):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("eta")
        case "E_ap":
            try:
                temp = float(input("Set the applied potential (V):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("E_ap")
    return temp

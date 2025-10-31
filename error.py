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
        case "ko":
            try:
                temp = float(input("Set ko:\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("ko")
        case "eta":
            try:
                temp = float(input("Set the overpotential:\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("eta")
        case "E_ap":
            try:
                temp = float(input("Set the applied potential:\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("E_ap")
    return temp

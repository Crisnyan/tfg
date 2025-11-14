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
        case "n_step":
            try:
                temp = int(input("Set the number of steps calculated:\n"))
                while temp < 10:
                    print('Error: the number of steps cannot be less than 10') 
                    temp = error("n_step")
            except Exception as e:
                print("Error:", e)
                temp = error("n_step")
        case "temperature":
            try:
                temp = float(input("Set the desired temperature (K):\n"))
                if temp == 0:
                    temp = 1e-12
                while temp < 0:
                    print('Error: the temperature cannot be negative') 
                    temp = error("temperature")
            except Exception as e:
                print("Error:", e)
                temp = error("temperature")
        case "endtime":
            try:
                temp = int(input("Set time limit in seconds (0 for no time limit):\n"))
                while temp < 0:
                    print('Error: time cannot be negative') 
                    temp = error("endtime")
            except Exception as e:
                print("Error:", e)
                temp = error("endtime")
        case "n_cycles":
            try:
                temp = int(input("Set the desired number of cycles:\n"))
                while temp <= 0:
                    print('Error: the number of cycles should be at least 1') 
                    temp = error("n_cycles")
            except Exception as e:
                print("Error:", e)
                temp = error("n_cycles")
        case "srate":
            try:
                temp = float(input("Set scan rate (V/s):\n"))
                while temp < 0:
                    print('Error: scan rate cannot be negative') 
                    temp = error("srate")
            except Exception as e:
                print("Error:", e)
                temp = error("srate")
        case "diff":
            try:
                temp = float(input("Set the diffusivity (m^2/s):\n"))
                while temp < 0:
                    print('Error: the diffusivity cannot be negative') 
                    temp = error("diff")
            except Exception as e:
                print("Error:", e)
                temp = error("diff")
        case "Cob":
            try:
                temp = float(input("Set the concentration of the oxidized species:\n"))
                if temp == 0:
                    temp = 1.0
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("Cob")
            except Exception as e:
                print("Error:", e)
                temp = error("Cob")
        case "Crb":
            try:
                temp = float(input("Set the concentration of the reduced species:\n"))
                if temp == 0:
                    temp = 1.0
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("Crb")
            except Exception as e:
                print("Error:", e)
                temp = error("Crb")
        case "c_red_an":
            try:
                temp = float(input("Set the concentration of the reduced species in the anode (1 for solids):\n"))
                if temp == 0:
                    temp = 1.0
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
                    temp = 1.0
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
                    temp = 1.0
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
                    temp = 1.0
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("c_ox_cat")
            except Exception as e:
                print("Error:", e)
                temp = error("c_ox_cat")
        case "intensity":
            try:
                temp = float(input("Set the intensity of the battery (A):\n"))
                if temp == 0:
                    temp = 1e-12
                if temp < 0:
                    temp *= -1
            except Exception as e:
                print("Error:", e)
                temp = error("intensity")
        case "resistance":
            try:
                temp = float(input("Set the internal resistance of the battery (ohms):\n"))
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
                temp = float(input("Set the nominal capacity of the battery (A/h):\n"))
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
                temp = float(input("Set initial jo (A/m^2):\n"))
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
                    temp = error("delta")
                while temp > 1:
                    print('Error: The diffusion layer must be 1 or less') 
                    temp = error("delta")
            except Exception as e:
                print("Error:", e)
                temp = error("delta")
        case "SOC":
            try:
                temp = float(input("Set state of charge of the battery (%):\n"))
                while temp < 0:
                    print('Error: The state of charge of the battery cannot be negative') 
                    temp = error("SOC")
                while temp > 100:
                    print('Error: The state of charge must be 100% or less') 
                    temp = error("SOC")
            except Exception as e:
                print("Error:", e)
                temp = error("SOC")
        case "k0":
            try:
                temp = float(input("Set k0 (mol/m^3):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("k0")
        case "eta":
            try:
                temp = float(input("Set the overpotential (V):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("eta")
        case "E_start":
            try:
                temp = float(input("Set the starting potential for the sweep (V):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("E_start")
        case "E_vertex":
            try:
                temp = float(input("Set the peak potential for the sweep (V):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("E_vertex")
        case "E_ap":
            try:
                temp = float(input("Set the applied potential (V):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("E_ap")
        case "E_ref":
            try:
                temp = float(input("Set the potential of the reference electrode, 0V for SHE (V):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("E_ap")
    return temp

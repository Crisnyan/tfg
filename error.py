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
        case "rpm_max":
            try:
                temp = int(input("Set the max RPM:\n"))
                while temp < 500:
                    print('Error: the max RPM must be at least 500') 
                    temp = error("rpm_max")
            except Exception as e:
                print("Error:", e)
                temp = error("rpm_max")
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
        case "viscosity":
            try:
                temp = float(input("Set the viscosity (m^2/s):\n"))
                while temp < 0:
                    print('Error: the viscosity cannot be negative') 
                    temp = error("viscosity")
            except Exception as e:
                print("Error:", e)
                temp = error("viscosity")
        case "area":
            try:
                temp = float(input("Set the electrode area (m^2):\n"))
                if temp == 0:
                    print('Error: area cannot be 0, set to 1') 
                    temp = 1.0
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("area")
            except Exception as e:
                print("Error:", e)
                temp = error("area")
        case "C0":
            try:
                temp = float(input("Set the initial concentration (mol/m^3):\n"))
                if temp == 0:
                    temp = 1.0
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("C0")
            except Exception as e:
                print("Error:", e)
                temp = error("C0")
        case "Cob":
            try:
                temp = float(input("Set the concentration of the oxidized species (mol/m^3):\n"))
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
                temp = float(input("Set the concentration of the reduced species (mol/m^3):\n"))
                if temp == 0:
                    temp = 1.0
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("Crb")
            except Exception as e:
                print("Error:", e)
                temp = error("Crb")
        case "c_red":
            try:
                temp = float(input("Set the concentration of the reduced species in mol/m^3 (1 for solids):\n"))
                if temp == 0:
                    temp = 1.0
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("c_red")
            except Exception as e:
                print("Error:", e)
                temp = error("c_red")
        case "c_ox":
            try:
                temp = float(input("Set the concentration of the oxidized species in mol/m^3 (1 for solids):\n"))
                if temp == 0:
                    temp = 1.0
                while temp < 0:
                    print('Error: concentration cannot be negative') 
                    temp = error("c_ox")
            except Exception as e:
                print("Error:", e)
                temp = error("c_ox")
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
                temp = float(input("Set the nominal capacity of the battery (A/s):\n"))
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
        case "j0":
            try:
                temp = float(input("Set j0 (A/m^2):\n"))
            except Exception as e:
                print("Error:", e)
                temp = error("j0")
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
                temp = error("E_ref")
        case _:
            print(f"ERROR, unbound")
                
    return temp

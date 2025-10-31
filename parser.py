
stdRedPotFile : dict[str, float] = {}
BatteryFile : dict[str, float] = {}

def parse(name: str) -> dict[str, float]: 
    d: dict[str, float] = {}
    try:
        with open(name, "r") as f:
            next(f)
            for line in f:
                key, value = [s.strip() for s in line.split(",")]
                d[key] = float(value)
            return d

    except OSError:
        print(f'Could not read file: ', name)
        exit()


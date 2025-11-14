
stdRedPotFile : dict[str, float] = {}
BatteryValuesFile : dict[float, float] = {}

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

def convert(strdict: dict[str, float]) -> dict[float, float]:
    d: dict[float, float] = {}
    for key, value in strdict.items():
        d[float(key)] = value
    return d

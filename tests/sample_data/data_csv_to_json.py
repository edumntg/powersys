import json
import pandas as pd
import numpy as np
# load buses data
buses = pd.read_csv('ieee9_buses.csv').to_numpy()
lines = pd.read_csv('ieee9_lines.csv').to_numpy()
generators = pd.read_csv('ieee9_gens.csv').to_numpy()

buses_json = {}
lines_json = {}
generators_json = {}

for row in buses:
    buses_json[int(row[0])] = {
        "V": row[1],"angle": row[2],"Pg": row[3],
        "Qg": row[4],
        "Pl": row[5],
        "Ql": row[6],
        "Vmin": row[7],
        "Vmax": row[8]
    }

for row in lines:
    lines_json[int(row[0])] = {
        "from": int(row[1]),
        "to": int(row[2]),
        "R": row[3],
        "X": row[4],
        "B": row[5],
        "a": row[6],
        "max_mva": row[7]
    }

for row in generators:
    generators_json[int(row[0])] = {
        "bus": int(row[1]),
        "c": row[2],
        "b": row[3],
        "a": row[4],
        "Pmin": row[5],
        "Pmax": row[6],
        "Qmin": row[7],
        "Qmax": row[8],
        "Ra": row[9],
        "Xd": row[10],
        "Xdp": row[12],
        "Xdpp": row[14],
        "Xq": row[11],
        "Xqp": row[13],
        "Xqpp": row[15],
        "Td0p": row[16],
        "Tq0p": row[17],
        "Td0pp": row[18],
        "Tq0pp": row[19],
        "H": row[20],
        "USE_CAGE": int(row[21]),
        "USE_AGC": int(row[22]),
        "Kv": row[23],
        "Ki": row[24],
        "Kt": row[25],
        "R": row[26],
        "Tv": row[27],
        "Tt": row[28],
        "USE_AVR": int(row[29]),
        "Kmed": row[30],
        "Kexc": row[31],
        "Ka": row[32],
        "Tmed": row[33],
        "Texc": row[34],
        "Ta": row[35],
        "Kd": row[36],
        "Kp": row[37],
        "Kvi": row[38],
        "Vexc_min": row[39],
        "Vexc_max": row[40],
        "USE_PSS": int(row[41]),
        "Kest": row[42],
        "Tw": row[43],
        "T1": row[44],
        "T2": row[45],
        "T3": row[46],
        "T4": row[47],
        "Vpc2_min": row[48],
        "Vpc2_max": row[49]
    }

final_json = {
    "buses": buses_json,
    "lines": lines_json,
    "generators": generators_json
}

json_object = json.dumps(final_json, indent=4)
with open("ieee9_buses.json", "w") as file:
    file.write(json_object)
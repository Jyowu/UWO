#!/usr/bin/env python3
"""
Generate constant/boundaryData/inlet/ files for timeVaryingMappedFixedValue
from the Figure 7 radial velocity profile (Dellenback et al., S=0.6).

Run from the case directory after running:
    postProcess -func writeCellCentres -time 0
"""

import os
import re
import subprocess
import numpy as np

# ── Digitised from Figure 7 (S=0.6, Dellenback et al.) ───────────────────────
# Axial velocity profile (r/R, Uz [m/s])
r_axial = np.array([
    0.0040, 0.0245, 0.0425, 0.0630, 0.0809, 0.0988, 0.1194, 0.1399, 0.1617,
    0.1861, 0.2067, 0.2260, 0.2453, 0.2620, 0.2826, 0.3032, 0.3224, 0.3405,
    0.3585, 0.3752, 0.3958, 0.4163, 0.4407, 0.4639, 0.4883, 0.5127, 0.5358,
    0.5602, 0.5872, 0.6090, 0.6308, 0.6500, 0.6705, 0.7000, 0.7282, 0.7487,
    0.7704, 0.7897, 0.8114, 0.8338, 0.8574, 0.8714, 0.8861, 0.9021, 0.9129,
    0.9250, 0.9327, 0.9441, 0.9530, 0.9682, 0.9771, 0.9833, 0.9896, 0.9958,
    0.9982, 0.9994, 0.9992, 0.9978, 0.9976, 0.9987, 1.0012, 1.0010, 0.9996,
    1.0020, 1.0018,
])
Uz_profile = np.array([
    0.2304, 0.2304, 0.2304, 0.2329, 0.2329, 0.2378, 0.2439, 0.2526, 0.2637,
    0.2747, 0.2871, 0.2994, 0.3129, 0.3265, 0.3400, 0.3536, 0.3647, 0.3782,
    0.3930, 0.4041, 0.4164, 0.4300, 0.4423, 0.4559, 0.4706, 0.4860, 0.4940,
    0.5027, 0.5125, 0.5211, 0.5273, 0.5298, 0.5347, 0.5366, 0.5384, 0.5408,
    0.5409, 0.5433, 0.5378, 0.5316, 0.5248, 0.5187, 0.5125, 0.5070, 0.4984,
    0.4903, 0.4793, 0.4694, 0.4546, 0.4374, 0.4164, 0.3967, 0.3733, 0.3499,
    0.3253, 0.3031, 0.2772, 0.2489, 0.2267, 0.2033, 0.1786, 0.1503, 0.1318,
    0.1072, 0.0850,
])

# Tangential velocity profile (r/R, Utheta [m/s])
r_tang = np.array([
    0.0018, 0.0148, 0.0239, 0.0368, 0.0498, 0.0614, 0.0693, 0.0809, 0.0951,
    0.1043, 0.1210, 0.1352, 0.1520, 0.1688, 0.1804, 0.1933, 0.2152, 0.2320,
    0.2526, 0.2693, 0.2912, 0.3105, 0.3311, 0.3478, 0.3697, 0.3889, 0.4159,
    0.4390, 0.4659, 0.4890, 0.5184, 0.5427, 0.5619, 0.5836, 0.6092, 0.6335,
    0.6526, 0.6756, 0.7011, 0.7203, 0.7433, 0.7713, 0.7854, 0.8148, 0.8313,
    0.8581, 0.8772, 0.8989, 0.9180, 0.9397, 0.9537, 0.9728, 0.9829, 0.9891,
    0.9940, 0.9965, 0.9976, 0.9975, 1.0025, 0.9984, 1.0009, 1.0007, 1.0019,
])
Utheta_profile = np.array([
    0.0862, 0.1047, 0.1306, 0.1528, 0.1725, 0.1947, 0.2205, 0.2402, 0.2600,
    0.2846, 0.3043, 0.3228, 0.3437, 0.3634, 0.3819, 0.3967, 0.4164, 0.4324,
    0.4485, 0.4645, 0.4805, 0.4990, 0.5138, 0.5236, 0.5409, 0.5520, 0.5643,
    0.5704, 0.5766, 0.5803, 0.5803, 0.5754, 0.5692, 0.5643, 0.5544, 0.5446,
    0.5359, 0.5248, 0.5150, 0.5039, 0.4965, 0.4867, 0.4817, 0.4682, 0.4534,
    0.4349, 0.4201, 0.4041, 0.3893, 0.3745, 0.3622, 0.3425, 0.3228, 0.2969,
    0.2747, 0.2526, 0.2279, 0.2033, 0.1823, 0.1565, 0.1318, 0.1035, 0.0899,
])

# ── Geometry ──────────────────────────────────────────────────────────────────
R = 0.05078 / 2.0   # inlet pipe radius [m]

# ── Read face centres from writeCellCentres output ───────────────────────────
def read_boundary_values(filepath, patch="inlet"):
    result = subprocess.run(
        ['foamDictionary', filepath, '-entry', f'boundaryField/{patch}/value', '-value'],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise ValueError(f"foamDictionary failed on {filepath}:\n{result.stderr}")
    text = result.stdout
    m = re.search(r'nonuniform List<scalar>\s+(\d+)\s*\((.*?)\)', text, re.DOTALL)
    if not m:
        raise ValueError(f"Could not parse output for patch '{patch}' from {filepath}")
    n = int(m.group(1))
    vals = list(map(float, m.group(2).split()))
    assert len(vals) == n, f"Expected {n} values, got {len(vals)}"
    return np.array(vals)

print("Reading face centres...")
cx = read_boundary_values("0/Ccx", patch="inlet")
cy = read_boundary_values("0/Ccy", patch="inlet")
cz = read_boundary_values("0/Ccz", patch="inlet")
n_faces = len(cx)
print(f"  Found {n_faces} inlet faces")
print(f"  z range: {cz.min():.5f} to {cz.max():.5f}")

# ── Compute radial distance and interpolate profiles ──────────────────────────
r = np.sqrt(cx**2 + cy**2)
r_norm = r / R
r_norm = np.clip(r_norm, 0.0, 1.0)

Uz  = np.interp(r_norm, r_axial, Uz_profile)
Uth = np.interp(r_norm, r_tang,  Utheta_profile)

# Convert tangential to Cartesian (anti-clockwise positive when viewed from +z)
# Uθ acts perpendicular to r in the x-y plane: (-y/r, x/r, 0)
with np.errstate(invalid="ignore", divide="ignore"):
    Ux = np.where(r > 1e-10, -Uth * cy / r, 0.0)
    Uy = np.where(r > 1e-10,  Uth * cx / r, 0.0)

# ── Write boundaryData files ──────────────────────────────────────────────────
out_dir = os.path.join("constant", "boundaryData", "inlet")
time_dir = os.path.join(out_dir, "0")
os.makedirs(time_dir, exist_ok=True)

# points file — raw list, no FoamFile header (OF12 reads directly as stream)
points_path = os.path.join(out_dir, "points")
with open(points_path, "w") as f:
    f.write(f"{n_faces}\n(\n")
    for i in range(n_faces):
        f.write(f"({cx[i]:.8f} {cy[i]:.8f} {cz[i]:.8f})\n")
    f.write(")\n")
print(f"Written: {points_path}")

# U file — raw list, no FoamFile header
U_path = os.path.join(time_dir, "U")
with open(U_path, "w") as f:
    f.write(f"{n_faces}\n(\n")
    for i in range(n_faces):
        f.write(f"({Ux[i]:.6f} {Uy[i]:.6f} {Uz[i]:.6f})\n")
    f.write(")\n")
print(f"Written: {U_path}")

print("\nProfile statistics:")
print(f"  Uz  : min={Uz.min():.3f}  max={Uz.max():.3f}  mean={Uz.mean():.3f} m/s")
print(f"  Uth : min={Uth.min():.3f}  max={Uth.max():.3f}  mean={Uth.mean():.3f} m/s")
print("\nDone. Now update 0/U inlet BC to timeVaryingMappedFixedValue.")

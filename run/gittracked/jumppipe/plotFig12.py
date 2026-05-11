#!/usr/bin/env python3
"""
Plot Figure 12: radial profiles of axial and tangential velocity
at z/D = 0.5, 1, 3 downstream of the expansion (S=0.6).

Run after simulation completes:
    postProcess -func sample -latestTime        (serial)
    postProcess -func sample -latestTime -parallel  (parallel, before reconstructPar)

Then: python3 plotFig12.py
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Dellenback et al. (1988) experimental data ────────────────────────────────
exp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Experimental_Data')

def _load_exp(fname):
    return np.loadtxt(os.path.join(exp_dir, fname), delimiter=',')

exp_uz = {
    0.5: _load_exp('zbyD0.5_Axial.csv'),
    1.0: _load_exp('zbyD1_Axial.csv'),
    3.0: _load_exp('zbyD3_Axial.csv'),
}
exp_uth = {
    0.5: _load_exp('zbyD0.5_Tangential.csv'),
    1.0: _load_exp('zbyD1_Tangential.csv'),
    3.0: _load_exp('zbyD3_Tangential.csv'),
}

# ── Case directories ──────────────────────────────────────────────────────────
cases = {
    'kOmegaSST':    {'label': 'SST',     'color': 'steelblue',  'ls': '--', 'enabled': True},
    'kOmegaSSTCC':  {'label': 'SST-CC',  'color': 'darkorange', 'ls': '-.', 'enabled': True},
    'kOmegaSSTCCS': {'label': 'SST-CCS', 'color': 'green',      'ls': '-',  'enabled': True},
    'realKE':       {'label': 'RKE',     'color': 'red',        'ls': ':',  'enabled': True},
    'realKECC':     {'label': 'RKE-CC',  'color': 'purple',     'ls': '--', 'enabled': True},
    'realKECCS':    {'label': 'RKE-CCS', 'color': 'brown',      'ls': '-.', 'enabled': True},
    'RSM':          {'label': 'RSM',     'color': 'black',      'ls': '-',  'enabled': True},
}
cases = {k: v for k, v in cases.items() if v['enabled']}

base_dir = os.path.dirname(os.path.abspath(__file__))

z_labels  = [0.5, 1.0, 3.0]
set_names = ['r_zD0p5', 'r_zD1', 'r_zD3']

PLOT_TIME = None   # set to None to use latest available time

# ── Read OF sample data ───────────────────────────────────────────────────────
def read_sample(case_dir, set_name, plot_time=None):
    """Return (r, Ux, Uy, Uz) arrays from postProcessing/sample.
    File format columns: x  U_x  U_y  U_z  (header line starting with #)
    """
    if plot_time is not None:
        pattern = os.path.join(case_dir, 'postProcessing', 'sample',
                               str(plot_time), f'{set_name}.xy')
        files = glob.glob(pattern)
    else:
        pattern = os.path.join(case_dir, 'postProcessing', 'sample', '*',
                               f'{set_name}.xy')
        files = glob.glob(pattern)
        files.sort(key=lambda p: float(os.path.basename(os.path.dirname(p))))
    if not files:
        return None
    data = np.loadtxt(files[-1], comments='#')
    r, Ux, Uy, Uz = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    return r, Ux, Uy, Uz

# ── R² ────────────────────────────────────────────────────────────────────────
def compute_r2(r_cfd, u_cfd, r_exp, u_exp):
    u_interp = np.interp(r_exp, r_cfd, u_cfd)
    ss_res = np.sum((u_exp - u_interp) ** 2)
    ss_tot = np.sum((u_exp - np.mean(u_exp)) ** 2)
    return np.nan if ss_tot == 0 else 1.0 - ss_res / ss_tot

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 2, figsize=(10, 12))
fig.subplots_adjust(hspace=0.35, wspace=0.3)
r2_scores = {case: {'uz': {}, 'uth': {}} for case in cases}

for row, (zD, sname) in enumerate(zip(z_labels, set_names)):
    ax_uz  = axes[row, 0]
    ax_uth = axes[row, 1]

    # Experimental data
    if exp_uz[zD].size > 0:
        r_e, uz_e = exp_uz[zD][:, 0], exp_uz[zD][:, 1]
        ax_uz.plot(r_e, uz_e, 'o', color='k', ms=4, label='Exp.', zorder=5)
    if exp_uth[zD].size > 0:
        r_e, uth_e = exp_uth[zD][:, 0], exp_uth[zD][:, 1]
        ax_uth.plot(r_e, uth_e, 'o', color='k', ms=4, label='Exp.', zorder=5)

    R2 = 0.04926   # outer pipe radius [m]

    # CFD results
    for case, props in cases.items():
        case_dir = os.path.join(base_dir, case)
        result = read_sample(case_dir, sname, plot_time=PLOT_TIME)
        if result is None:
            continue
        x, Ux, Uy, Uz = result
        # Shift x so left wall=0, centreline=R2, right wall=2R2=D2
        x_plot = x + R2
        # Tangential: negated Uy to match Dellenback Fig 12 sign convention
        Uth = -Uy
        ax_uz.plot(x_plot, Uz, color=props['color'], ls=props['ls'],
                   lw=1.5, label=props['label'])
        ax_uth.plot(x_plot, Uth, color=props['color'], ls=props['ls'],
                    lw=1.5, label=props['label'])
        r2_scores[case]['uz'][zD]  = compute_r2(x_plot, Uz,  exp_uz[zD][:, 0],  exp_uz[zD][:, 1])
        r2_scores[case]['uth'][zD] = compute_r2(x_plot, Uth, exp_uth[zD][:, 0], exp_uth[zD][:, 1])

    ax_uz.set_title(f'z/D = {zD}', fontsize=10)
    ax_uth.set_title(f'z/D = {zD}', fontsize=10)
    ax_uz.set_ylabel('Axial velocity (m/s)')
    ax_uth.set_ylabel('Tangential velocity (m/s)')
    ax_uz.axhline(0, color='k', lw=0.5, ls=':')
    ax_uth.axhline(0, color='k', lw=0.5, ls=':')
    ax_uz.set_xlim(0, 2*R2)
    ax_uth.set_xlim(0, 2*R2)
    ax_uz.set_ylim(-0.3, 0.9)
    ax_uth.set_ylim(-0.5, 0.5)

for ax in axes[-1, :]:
    ax.set_xlabel('Radial distance (m)')

# Legend on bottom-right subplot
handles, labels = axes[1, 1].get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, loc='lower center', ncol=4,
               bbox_to_anchor=(0.5, -0.02), fontsize=8)

print("\nR² scores vs Dellenback et al. experimental data")
header = f"{'Case':<14}" + "".join(f"  z/D={zD:<3}Uz    Uth  " for zD in z_labels)
print(header)
for case, scores in r2_scores.items():
    row = f"{cases[case]['label']:<14}"
    for zD in z_labels:
        uz  = scores['uz'].get(zD, float('nan'))
        uth = scores['uth'].get(zD, float('nan'))
        row += f"  {uz:6.3f}  {uth:6.3f}  "
    print(row)

plt.suptitle('Mean velocity profiles — S=0.6 swirling flow (Dellenback et al.)',
             fontsize=11)
out = f'Velocity_profiles_t{PLOT_TIME if PLOT_TIME is not None else "latest"}.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f'Saved: {out}')
plt.close()

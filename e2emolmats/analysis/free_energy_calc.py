"""
To calculate free energy difference between acridine melt and crystal phases.

\Delta G(T) = (T_m-T) * (\Delta H_m/T_m - (C_(p,l0) - C_(p,s0)) + (T_m+3*T)/2 + T*(C_(p,l0) - C_(p,s0) * ln(T_m/T))
"""

import numpy as np
import plotly.graph_objects as go

# experimental data
polymorphs_dict = {
    'Form 3': {
        'T_melt': 108.2,  # K
        'H_fus': 20.41,  # kJ/mol
        'C_p1': 0.205,  # kJ/mol*K
        'C_p0': 0  # kJ/mol*K
    },
    'Form 4': {
        'T_melt': 110.1,  # K
        'H_fus': 21.01,  # kJ/mol
        'C_p1': 0.205,  # kJ/mol*K
        'C_p0': 0  # kJ/mol*K
    },
    'Form 6': {
        'T_melt': 108.5,  # K
        'H_fus': 19.69,  # kJ/mol
        'C_p1': 0.205,  # kJ/mol*K
        'C_p0': 0  # kJ/mol*K
    },
    'Form 7': {
        'T_melt': 107.6,  # K
        'H_fus': 20.41,  # kJ/mol
        'C_p1': 0.205,  # kJ/mol*K
        'C_p0': 0  # kJ/mol*K
    },
    'Form 8': {
        'T_melt': 108.8,  # K
        'H_fus': 22.03,  # kJ/mol
        'C_p1': 0.205,  # kJ/mol*K
        'C_p0': 0  # kJ/mol*K
    },
}

num_polymorphs = len(polymorphs_dict)
t_range = np.linspace(273, 420, 100)

deltaG = np.zeros((len(t_range), num_polymorphs))

cpl1 = 0.0205
cpl0 = 0

for ind, T in enumerate(t_range):
    for ind2, (polymorph, values) in enumerate(polymorphs_dict.items()):
        cp1 = values['C_p1']
        cp0 = values['C_p0']
        Tmelt = values['T_melt']
        latent = values['H_fus']

        deltaG[ind, ind2] = (T - Tmelt) * (
                latent / Tmelt - (cpl0 - cp0) + ((Tmelt + 3 * T) / 2) * (cpl1 - cp1) + T * (cpl0 - cp0) * np.log(
            Tmelt / T)
        )

fig = go.Figure()
for ind, polymorph in enumerate(polymorphs_dict.keys()):
    fig.add_scattergl(x=t_range, y=deltaG[:, 1] - deltaG[:, ind], name=polymorph, mode='lines')

fig.show(renderer='browser')

aa = 1

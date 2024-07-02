"""
To calculate free energy difference between acridine melt and crystal phases.

\Delta G(T) = (T_m-T) * (\Delta H_m/T_m - (C_(p,l0) - C_(p,s0)) + (T_m+3*T)/2 + T*(C_(p,l0) - C_(p,s0) * ln(T_m/T))
"""

import numpy as np
import plotly.graph_objects as go

# experimental data
experimental_polymorphs_dict = {
    'Form3': {
        'T_melt': 108.2 + 273.15,  # K
        'H_fus': 20.41,  # kJ/mol
    },
    'Form4': {
        'T_melt': 110.1 + 273.15,  # K
        'H_fus': 21.01,  # kJ/mol
    },
    'Form6': {
        'T_melt': 108.5 + 273.15,  # K
        'H_fus': 19.69,  # kJ/mol
        'C_p1': 0.205,  # kJ/mol*K
        'C_p0': 0  # kJ/mol*K
    },
    'Form7': {
        'T_melt': 107.6 + 273.15,  # K
        'H_fus': 20.41,  # kJ/mol
    },
    'Form8': {
        'T_melt': 108.8 + 273.15,  # K
        'H_fus': 22.03,  # kJ/mol
    },
}

cpl0_min_cps0 = 0.14463
cpl1_min_cps1 = -0.000272

theoretical_polymorphs_dict = {
    'Form2': {
        'T_melt': 395.23,

    },
    'Form3': {
        'T_melt': 394.68,  # K

    },
    'Form4': {
        'T_melt': 351.35,  # K

    },
    'Form6': {
        'T_melt': 358.82,  # K

    },
    'Form7': {
        'T_melt': 372.52,  # K

    },
    'Form8': {
        'T_melt': 371.71,  # K

    },
    'Form9': {
        'T_melt': 376.03,

    },
}
calculated_cp_coeffs = {'Form2': ([-1.58729441e-04,  6.48404836e-01]),
 'Form3': ([-7.44433207e-05,  6.24923029e-01]),
 'Form4': ([-0.00152371,  1.10869871]),
 'Form6': ([0.00174128, 0.01547051]),
 'Form7': ([0.00096636, 0.25308288]),
 'Form8': ([-3.99229567e-04,  7.27688252e-01]),
 'Form9': ([ 0.00208767, -0.15329192])}

calculated_latents = {'Form2': 21.235253259744866,
 'Form3': 18.44023695338433,
 'Form4': 14.550585862221197,
 'Form6': 15.2178183848456,
 'Form7': 18.484026829707744,
 'Form8': 15.393464084682563,
 'Form9': 19.558361932176666}

calculated_melt_cp_coefficients = [-1.57578735e-05,  6.51564825e-01]

for key in theoretical_polymorphs_dict.keys():
    theoretical_polymorphs_dict[key]['H_fus'] = calculated_latents[key]
    theoretical_polymorphs_dict[key]['C_p1'] = calculated_cp_coeffs[key][0]
    theoretical_polymorphs_dict[key]['C_p0'] = calculated_cp_coeffs[key][1]

num_polymorphs = len(experimental_polymorphs_dict)
t_range = np.linspace(273, 420, 100)

deltaG_exp = np.zeros((len(t_range), num_polymorphs))
'''experimental values'''
for ind, temp in enumerate(t_range):
    for ind2, (polymorph, values) in enumerate(experimental_polymorphs_dict.items()):
        Tmelt = values['T_melt']
        latent = values['H_fus']

        deltaG_exp[ind, ind2] = (temp - Tmelt) * (
                (latent / Tmelt) - cpl0_min_cps0 + (
                    (Tmelt + 3 * temp) / 2) * cpl1_min_cps1) + temp * cpl0_min_cps0 * np.log(Tmelt / temp)

deltaG_comp = np.zeros((len(t_range), num_polymorphs))
cpl1 = calculated_melt_cp_coefficients[0]
cpl0 = calculated_melt_cp_coefficients[1]
'''computed values'''
for ind, temp in enumerate(t_range):
    for ind2, (polymorph, values) in enumerate(experimental_polymorphs_dict.items()):
        values = theoretical_polymorphs_dict[polymorph]
        cps1 = values['C_p1']
        cps0 = values['C_p0']
        Tmelt = values['T_melt']
        latent = values['H_fus']

        deltaG_comp[ind, ind2] = (temp - Tmelt) * (
                (latent / Tmelt) - (cpl0 - cps0) + (
                (Tmelt + 3 * temp) / 2) * (cpl1 - cps1)) + temp * (cpl0 - cps0) * np.log(Tmelt / temp)

from plotly.subplots import make_subplots
from plotly.colors import n_colors

colors = n_colors('rgb(25, 200, 255)', 'rgb(255, 0, 0)', len(experimental_polymorphs_dict.keys()), colortype='rgb')
fig = make_subplots(1, 2, subplot_titles=['Experimental', 'Simulations'])
for ind, polymorph in enumerate(experimental_polymorphs_dict.keys()):
    tmelt_ind = np.argmin(np.abs(experimental_polymorphs_dict[polymorph]['T_melt'] - t_range))
    fig.add_scattergl(x=t_range[:tmelt_ind] - 273.15, y=-deltaG_exp[:tmelt_ind, ind] + deltaG_exp[:tmelt_ind, 1],
                      name=polymorph, mode='lines',
                      marker_color=colors[ind],
                      row=1, col=1)

    tmelt_ind = np.argmin(np.abs(theoretical_polymorphs_dict[polymorph]['T_melt'] - t_range))
    fig.add_scattergl(x=t_range[:tmelt_ind] - 273.15, y=-deltaG_comp[:tmelt_ind, ind] + deltaG_comp[:tmelt_ind, 1],
                      name=polymorph, mode='lines',
                      marker_color=colors[ind],
                      showlegend=False,
                      row=1, col=2)

fig.update_layout(xaxis_title='Temperature (C)', yaxis_title='Free Energy (kJ/mol)')
fig.show(renderer='browser')

aa = 1

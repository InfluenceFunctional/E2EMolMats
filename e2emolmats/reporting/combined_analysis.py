'multi-battery analysis'
from e2emolmats.analysis.utils import atoms_per_molecule
from e2emolmats.reporting.temp_analysis import lattice_energy_figs, mean_temp_anomaly_fig, temperature_profile_fig, \
    com_deviation_fig, single_run_thermo_fig, ramped_melt_T_extraction
from e2emolmats.reporting.utils import cp_and_latent_analysis, confirm_melt, analyze_heat_capacity, \
    latent_heat_analysis, crystal_stability_analysis, compute_and_plot_melt_slopes, compute_and_plot_melt_slopes_com, \
    plot_melt_points, POLYMORPH_MELT_POINTS

import numpy as np


def combined_trajectory_analysis(config, combined_df, wandb):
    if config.compute_melt_temps:
        #
        # import plotly.graph_objects as go
        # from scipy.ndimage import gaussian_filter1d
        #
        # fig = go.Figure()
        #
        # good_inds = np.argwhere((combined_df['structure_identifier'] == 'acridine/Form2') *
        #                         (combined_df['temperature'] == 350)).flatten()
        #
        # for indyind, ind in enumerate(good_inds):
        #     fig.add_trace(
        #         go.Scattergl(x=combined_df['time step'][ind] / 1e6,
        #                      y=gaussian_filter1d(combined_df['E_pair'][ind] / combined_df['num_molecules'][ind], 5),
        #                      marker_color='blue', name='Form2',
        #                      showlegend=indyind == 0))
        #
        # good_inds = np.argwhere((combined_df['structure_identifier'] == 'acridine/Form4') *
        #                         (combined_df['temperature'] == 350)).flatten()
        #
        # for indyind, ind in enumerate(good_inds):
        #     fig.add_trace(
        #         go.Scattergl(x=combined_df['time step'][ind] / 1e6,
        #                      y=gaussian_filter1d(combined_df['E_pair'][ind] / combined_df['num_molecules'][ind], 5),
        #                      marker_color='red', name='Form4',
        #                      showlegend=indyind == 0))
        #
        # fig.show(renderer='browser')
        combined_df = confirm_melt(combined_df)  # TODO note for run 5 the equil time was hardcoded!
        print(f"{np.sum(combined_df['Melt Succeeded'] != True)} failed to melt!")
        combined_df.drop(index=np.argwhere(combined_df['Melt Succeeded'] != True).flatten(), inplace=True)
        combined_df.reset_index(drop=True, inplace=True)

        # temperature directional profile
        dev_slopes = []
        if True:  # False:
            mean_temp_anomaly_fig(combined_df)

            for r_ind in range(len(combined_df)):
                temperature_profile_fig(combined_df, r_ind, sigma_x=1, sigma_y=2, show_fig=True)
                fig, deviation, com_dev_slope = com_deviation_fig(combined_df, r_ind, show_fig=False)
                dev_slopes.append(com_dev_slope)

        combined_df['com_deviation_slope'] = dev_slopes

        fig, melt_estimate_dict, _ = compute_and_plot_melt_slopes(combined_df)
        fig, melt_estimate_dict2 = compute_and_plot_melt_slopes_com(combined_df)

        fig2 = plot_melt_points(melt_estimate_dict, POLYMORPH_MELT_POINTS[config.molecule])
        fig3 = plot_melt_points(melt_estimate_dict2, POLYMORPH_MELT_POINTS[config.molecule])

    if config.melt_scan_analysis:
        """
        We want several clean thermo analyses
    
        local T1-T3, pe, mobility
                all vs time (nicely break out simulation phases)
                all vs temperature
                all between bulk, crystal, interface
                
        then, we need to extract melt point as the temperature where melting starts
        """
        melt_temp_dict = {}
        for r_ind in range(len(combined_df)):
            row = combined_df.iloc[r_ind]
            single_run_thermo_fig(row)
            melting_temp = ramped_melt_T_extraction(row)
            melt_temp_dict[row['structure_identifier']] = melting_temp

        aa = 1

    if config.nanocluster_analysis:
        combined_df = confirm_melt(combined_df)
        combined_df.drop(index=np.argwhere(combined_df['Melt Succeeded'] != True).flatten(), inplace=True)
        combined_df.reset_index(drop=True, inplace=True)
        crystal_stability_analysis(combined_df)

    if config.latents_analysis:
        fig = latent_heat_analysis(combined_df)

    if config.cp_analysis:
        combined_df = analyze_heat_capacity(combined_df, atoms_per_molecule)

    if config.cp2_analysis:
        combined_df = confirm_melt(combined_df)
        cp_and_latent_analysis(combined_df).show(renderer='browser')

    if config.lattice_energy_analysis:
        lattice_energy_figs(combined_df)

    aa = 1



import os
import pandas as pd
import numpy as np

from reporting.cluster_figs import cluster_property_heatmap
from utils import compute_rdf_distance

os.chdir(r'C:\Users\mikem\crystals\clusters\cluster_structures')

# results_df = pd.read_pickle('pure_clusters_1/results_df')
# results_df = pd.concat([results_df, pd.read_pickle('defect_clusters_1/results_df')])
# results_df = pd.concat([results_df, pd.read_pickle('defect_clusters_2/results_df')])
#
# '''alignment analysis'''
# results_df = pd.read_pickle('defect_clusters_1/results_df')
# results_df = pd.concat([results_df, pd.read_pickle('defect_clusters_2/results_df')])
#
# cluster_property_heatmap(results_df, 'global_Ip_alignment', 'temperature', 'defect_rate',
#                          extra_axes=['max_sphere_radius'], extra_axes_values=[10], take_mean=True)
#
# cluster_property_heatmap(results_df, 'global_Ip_alignment', 'temperature', 'defect_rate',
#                          extra_axes=['max_sphere_radius'], extra_axes_values=[15], take_mean=True)
#
# cluster_property_heatmap(results_df, 'local_Ip_alignment', 'temperature', 'defect_rate',
#                          extra_axes=['max_sphere_radius'], extra_axes_values=[10], take_mean=True)
#
# cluster_property_heatmap(results_df, 'local_Ip_alignment', 'temperature', 'defect_rate',
#                          extra_axes=['max_sphere_radius'], extra_axes_values=[15], take_mean=True)
#
# del results_df
'''rdf analysis'''
results_df = pd.read_pickle('defect_clusters_1/results_df3')
results_df = pd.concat([results_df, pd.read_pickle('defect_clusters_2/results_df3')])
reference_df = pd.read_pickle('bulk_reference/reference_df3')
results_df = results_df.reset_index()
reference_df = reference_df.reset_index()

sample_ref_dists = np.zeros(len(results_df))
for i in range(len(results_df)):
    for j in range(len(reference_df)):
        if reference_df['temperature'][j] == 200 \
                and reference_df['structure_identifier'][j] == results_df['structure_identifier'][i]:
            ref_ind = j
            pass

    sample_rdf = results_df['intermolecular_rdfs'].loc[i][5:]
    reference_rdf = reference_df['intermolecular_rdfs'].loc[ref_ind][5:] # average over reference trajectory
    distmat = np.zeros((len(sample_rdf), len(reference_rdf)))
    for t1 in range(len(sample_rdf)):
        for t2 in range(len(reference_rdf)):
            distmat[t1,t2] = compute_rdf_distance(sample_rdf[t1], reference_rdf[t2], np.linspace(0,5,100), envelope='tanh')

    sample_ref_dists[i] = distmat.mean()

results_df['rdf_dists'] = sample_ref_dists

cluster_property_heatmap(results_df, 'rdf_dists', 'temperature', 'defect_rate',
                         extra_axes=['max_sphere_radius'], extra_axes_values=[10], take_mean=False)

cluster_property_heatmap(results_df, 'rdf_dists', 'temperature', 'defect_rate',
                         extra_axes=['max_sphere_radius'], extra_axes_values=[15], take_mean=False)

aa = 1
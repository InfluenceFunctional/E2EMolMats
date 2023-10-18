import os
import pandas as pd
import numpy as np

from reporting.cluster_figs import cluster_property_heatmap, collate_property_over_multiple_runs, plot_classifier_pies
from utils import compute_rdf_distance
from plotly.subplots import make_subplots
import plotly.graph_objects as go

os.chdir(r'C:\Users\mikem\crystals\clusters\cluster_structures')

results_df = pd.read_pickle('D:\crystals_extra\defect_clusters_5/results_df')
reference_df = pd.read_pickle('bulk_reference/reference_df5')
results_df = results_df.reset_index()
reference_df = reference_df.reset_index()

'''rdf analysis'''
sample_ref_dists = np.zeros(len(results_df))
for i in range(len(results_df)):
    for j in range(len(reference_df)):
        if reference_df['temperature'][j] == 200 \
                and reference_df['structure_identifier'][j] == results_df['structure_identifier'][i]:
            ref_ind = j
            pass

    sample_rdf = results_df['intermolecular_rdfs'].loc[i][5:]
    reference_rdf = reference_df['intermolecular_rdfs'].loc[ref_ind][5:]  # average over reference trajectory
    distmat = np.zeros((len(sample_rdf), len(reference_rdf)))
    for t1 in range(len(sample_rdf)):
        for t2 in range(len(reference_rdf)):
            distmat[t1, t2] = compute_rdf_distance(sample_rdf[t1], reference_rdf[t2], np.linspace(0, 5, 100), envelope='tanh')

    sample_ref_dists[i] = distmat.mean()

results_df['rdf_dists'] = sample_ref_dists

"""NN output analysis"""
classes = results_df['NN_classes'][0]
NNout_means = np.zeros((len(results_df), 10))
for i in range(len(results_df)):
    for j in range(len(classes)):
        # trailing 100 frames
        NNout_means[i, j] = np.mean(results_df['NN_trajectories'][i][-100:, j])

        # full trajectory mean
        # NNout_means[i, j] = results_df['NN_means'][i][j]

for i, label in enumerate(results_df['NN_classes'][0]):
    results_df[label] = NNout_means[:, i]

"""classifier outputs"""
plot_classifier_pies(results_df, 'max_sphere_radius', 'defect_rate', extra_axes=['temperature'], extra_axes_values=[250])

"alignment"
cluster_property_heatmap(results_df, 'global_Ip_alignment', 'max_sphere_radius', 'defect_rate',
                         extra_axes=['temperature'], extra_axes_values=[250], take_mean=True, norm_against_zero_y=False)

cluster_property_heatmap(results_df, 'local_Ip_alignment', 'max_sphere_radius', 'defect_rate',
                         extra_axes=['temperature'], extra_axes_values=[250], take_mean=True, norm_against_zero_y=False)

"RDF"
cluster_property_heatmap(results_df, 'rdf_dists',
                         'max_sphere_radius',
                         'defect_rate',
                         extra_axes=['temperature'], extra_axes_values=[250],
                         take_mean=False, norm_against_zero_y=False)

aa = 0

"""NN outputs"""
# ind = 0
#
# fig = make_subplots(rows=1, cols=3, subplot_titles=['Trajectory Mean', 'Variability', 'Trajectory'])
# categories = results_df['NN_classes'][ind]
# avg_probs = results_df['NN_means'][ind]
# atomwise_variance = results_df['NN_variance'][ind]
# prob_trajectories = results_df['NN_trajectories'][ind]
# fig.add_trace(go.Bar(x=categories, y=avg_probs, showlegend=False),
#               row=1, col=1)
# fig.add_trace(go.Bar(x=categories, y=atomwise_variance, showlegend=False),
#               row=1, col=2)
# for i in range(len(categories)):
#     fig.add_trace(go.Scattergl(y=prob_trajectories[:, i], name=categories[i]),
#                   row=1, col=3)
# fig.show(renderer='browser')

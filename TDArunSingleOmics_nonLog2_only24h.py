"""
Created on Thu Mar  2 10:09:05 2023

@author: Doroteya K. Staykova, Multicore Dynamics Ltd
"""
import os
import csv
import numpy as np
import pandas as pd
# TDA
from gtda.mapper import (
    CubicalCover,
    OneDimensionalCover,
    make_mapper_pipeline,
    plot_static_mapper_graph,
    plot_interactive_mapper_graph
)
from gtda.pipeline import make_pipeline
# Anndata
import anndata as ad
import scanpy as sc
# processing
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.cluster import DBSCAN
from timeit import default_timer as timer


class CompositeFilterFunc(BaseEstimator, TransformerMixin):
    def __init__(self, ann=[]):
        self.ann = ann
        pass

    def fit(self, X, y=None):
        self._is_fitted = True
        return self

    def transform(self, X, y=None):
        # %% PCoA
        embedding = self.ann.obsm['X_pcoa']

        lensFoC = self.ann.obs['Only last FoC measurment for each sample'].values.astype(
            float).reshape(-1, 1)
        lenses = np.c_[embedding, lensFoC]

        return lenses


def make_sample_name(idx):
    # Example:
    # WTC_B2_0h_C1_TR2
    # WTC_B2_24h_C3_TR1
    parts = idx.split("_")

    timepoint = next(p for p in parts if p.endswith("h"))  # 0h or 24h
    condition = next(p for p in parts if p.startswith("C"))
    replicate = next(p for p in parts if p.startswith("TR"))

    start, treatment = mapping[(timepoint, condition)]

    br = int(replicate.replace("TR", ""))
    eht = start + br - 1

    return f"EHT{eht}_{timepoint}_{treatment}_BR{br}"


# %%
# ====
# MAIN
# ====
randomValue = 42
np.random.seed(randomValue)
# %% load annotated data for 24h
flData = '24h_AnnData_20250222__Exp_EHTs-drugs__240808fc.h5ad'
dataPath = os.path.join(os.getcwd(), 'DATA', flData)
adata = sc.read_h5ad(dataPath)
dataTag = flData.split('.h5ad')[0]
# %% add updated labels (Fede: Supplementary data)
mapping = {
    ("0h", "C1"): (1, "UNTREAT"),
    ("24h", "C1"): (4, "UNTREAT"),
    ("24h", "C2"): (7, "DMSO"),
    ("24h", "C3"): (10, "DOXO"),
    ("24h", "C4"): (13, "EPI"),
}

adata.obs["Sample"] = adata.obs.index.map(make_sample_name)
# %% filter function(s)
col_lens = []

for i in range(2):
    lens_label = 'lens_PCoA%i' % (i+1)
    col_lens += [lens_label]

col_lens += ['lens_Only last FoC measurment for each sample']
netTag = '_PCoAplusFoC_i9frac08balanced'

inputData = adata.X
filter_func = CompositeFilterFunc(adata)
# %% TDA
startTime = timer()
'''Choose type of cover'''
_, projDim = inputData.shape

cover = CubicalCover(kind='balanced', n_intervals=9, overlap_frac=0.8)
'''Choose clustering algorithm'''
clusterer = DBSCAN(
    eps=0.5, min_samples=1, metric='correlation')  # metric applied to the original dataset
'''Setup TDA pipeline, the structure is adata → filter_func.transform → lenses → cover → clustering → graph
    1. filter_func (custom transformer) produces the lenses (low-dimensional representation)
    2. cover splits the lens space into overlapping regions (bins / hypercubes / balls)
    3. clusterer clusters data within each cover region
    4. Mapper builds the graph from overlapping clusters
'''
pipe = make_mapper_pipeline(
    filter_func=filter_func,
    cover=cover,
    clusterer=clusterer,
    verbose=True,
    n_jobs=-1,
)
'''for nodes enrichment with a categorical parameter'''
# combine 0h and non-treated
idx = adata.obs.index.tolist()
vlsNonTreated = list(map(lambda x: '%s_%s' % (
    adata[x].obs['Timepoint'].values[0], adata[x].obs['Condition'].values[0]) if adata[x].obs['Condition'].values[0] == 'Not treated' else 'Treated', idx))
adata.obs['State'] = vlsNonTreated

parOI = ['Condition']  # , 'State']
color_data_categorical = pd.get_dummies(
    adata.obs[parOI].astype('str'), prefix=parOI)
'''add numerical parameters'''
numerical_parOI = ['Relative FoC change']
color_data_num = pd.DataFrame()
for col in numerical_parOI:
    color_data_num[col] = adata.obs[col].astype(float)
color_data_num.index = color_data_categorical.index
'''add sample IDs'''
parOI = ['Sample']
color_data_id = pd.get_dummies(
    adata.obs[parOI].astype('str'), prefix=parOI)
color_data_id.index = color_data_num.index
'''add lenses'''
lenses = filter_func.fit_transform(inputData)
color_lens = pd.DataFrame(lenses, columns=col_lens)
color_lens.index = color_data_num.index
'''put together all data for colour'''
color_data = pd.concat(
    [color_lens, color_data_num, color_data_categorical, color_data_id], axis=1)
# %% colour choice
plot_choice = "inferno"
template_choice = "plotly_white"
plotly_params = {"node_trace": {"marker_colorscale": plot_choice},
                 "layout": {"template":  template_choice,
                            }}

fig = plot_static_mapper_graph(
    pipe, inputData, color_data=color_data, plotly_params=plotly_params)
# %% update layout
if template_choice == "plotly_white":
    fig.update_layout(plot_bgcolor="whitesmoke")
# %% display figure
fig.show(config={'scrollZoom': True})

procTime = timer() - startTime
print('\nTDA processing time: %.2f min' % round(procTime/60., 2))
# %%
exportHTML = True

if exportHTML:
    # ## Export the interactive network model as a standalone `*.html` file
    html_path = os.path.join(
        os.getcwd(), 'RESULTS', f"TCN_{netTag}.html")

    fig.write_html(
        html_path,
        full_html=True,
        auto_open=True,
        # config={"responsive": True},
    )
    print(f"Saved TDA network to:\n{html_path}")
# %%
mapper_graph = pipe.fit_transform(inputData)
print(type(mapper_graph))
print(mapper_graph.vs.attributes())

set_labels = mapper_graph.vs["pullback_set_label"]
node_elements = mapper_graph.vs["node_elements"]
tupleNodes = list(map(lambda x, y: (x, y), set_labels, node_elements))
sampleNames = adata.obs.index  # adata.obs['Sample']
# for each sample -> find associated nodes and save info
saveAllNodes = True
if saveAllNodes:
    nodeOutput = []
    uniqueCount = 0
    nodeHdr = ['unique ID', 'pullback_set_label',
               'node_elements', 'node_size']
    for elem in tupleNodes:
        nodeID = elem[0]
        uniqueCount += 1

        for k in elem[1]:  # each sample in node
            name = sampleNames[k]
            print('Node %i: s%s' % (nodeID, name))
            vecOut = [uniqueCount, nodeID, name, len(elem[1])]
            nodeOutput.append(vecOut)

    df_output = pd.DataFrame(nodeOutput, columns=nodeHdr)

    flOut = os.path.join(os.getcwd(), 'RESULTS',
                         'Nodes%s_%s.csv' % (netTag, dataTag))
    df_output.to_csv(flOut, index=False)

print('\nDone.')

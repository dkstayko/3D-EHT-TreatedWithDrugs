"""
Created on Fri Mar  7 10:19:55 2025

@author: Doroteya K. Staykova, Multicore Dynamics Ltd
"""
import os
import csv
import sys
import numpy as np
import pandas as pd
import scipy.stats as stats
# Anndata
import anndata as ad
import scanpy as sc
# timing
from timeit import default_timer as timer


def corr_fun(calc_type, meta_vls, vls):
    '''function to calculate correlation values
       args: calc_type - 'Pearson' or 'Spearman'
             meta_vls - metadata values to be correlated
             vls - numerical data coming from a column in dataframe (feature values)
    '''
    arr_vls = np.array(vls)
    arr_meta = np.array(meta_vls)
    # check if all feature values are 0
    if any([not np.any(arr_vls),
            not np.any(arr_meta),
            ]):
        return (0, 1)

    if calc_type == 'Pearson':
        PC_stat = stats.pearsonr(arr_meta, arr_vls)
        pc = PC_stat.statistic
        pc_p = PC_stat.pvalue

        return (pc, pc_p)

    if calc_type == 'Spearman':
        if not np.any(arr_vls):
            return (0, 0)
        '''Spearman correlation'''
        spearman_corr = stats.spearmanr(arr_meta, arr_vls)
        sc = spearman_corr.correlation
        sc_p = spearman_corr.pvalue

        return (sc, sc_p)


def cat2num(cat_val, val2compare):
    '''function to transform the categorical information into node enrichment value
        args: cat_val - categorical string information about node enrichment with ALL conditions
              val2compare - SINGLE condition of interest to calculate numerical node enrichment
    '''

    if val2compare not in cat_val:
        return 0.

    else:
        elems = cat_val.split('; ')
        for x in elems:
            if val2compare in x:
                enrichment = float(x.split(':')[-1])

                return enrichment


def get_single_protein_acc(vls):
    '''function to select the first of multiple values listed in protein_accession
       args: vls - protein_accession as coming from original data
    '''
    pa = vls.split('.')[1]
    pa_first = pa.split(';')[0]
    '''if dash found, remove as Reactome does not recognise isoforms'''
    if '-' in pa_first:
        # print('\nModified %s to %s' % (pa_first, pa_first.split('-')[0]))
        return pa_first.split('-')[0]

    else:
        return pa_first


# %%
# ====
# MAIN
# ====
randomValue = 42
np.random.seed(randomValue)
# %% load data and metadata
startTime = timer()

dataTag = '24h_AnnData_20250222__Exp_EHTs-drugs__240808fc'
dataFl = '%s.h5ad' % dataTag
dataPath = os.path.join(os.getcwd(), 'DATA', dataFl)
adata = sc.read_h5ad(dataPath)
print('\nSource data loaded.')
# %% calculate correlation of FoC with raw data
performRawCorr = True

if performRawCorr:
    meta = 'Only last FoC measurment for each sample'
    meta_vls = adata.obs[meta].values.astype(float)
    df = pd.DataFrame(adata.X, columns=adata.var_names)
    '''Pearson'''
    corrPearson = df.apply(lambda x: corr_fun(
        calc_type='Pearson', meta_vls=meta_vls, vls=x), axis=0)
    adata.varm['Pearson: %s' % meta] = corrPearson.iloc[0].values
    adata.varm['Pearson p-value: %s' % meta] = corrPearson.iloc[1].values
    '''Spearman'''
    corrSpearman = df.apply(lambda x: corr_fun(
        calc_type='Spearman', meta_vls=meta_vls, vls=x), axis=0)
    adata.varm['Spearman: %s' % meta] = corrSpearman.iloc[0].values
    adata.varm['Spearman p-value: %s' % meta] = corrSpearman.iloc[1].values
    ''''saving correlations as feature metadata in the trimmed data'''
    fl = os.path.join(os.getcwd(), 'RESULTS', 'Corr_%s.h5ad' % dataTag)
    adata.write(filename=fl)
    print('\nAnnotated data with correlations are saved in %s' % fl)

# %% load TCN info
nodeFl = 'Nodes_PCoAplusFoC_i9frac08balanced_%s.csv' % dataTag
tdaFl = 'TCN_PCoAplusFoC_i9frac08balanced.html'

rootPath = os.path.join(os.getcwd(), 'RESULTS')

dataPath = os.path.join(rootPath, nodeFl)
df_nodes = pd.read_csv(dataPath)
print('\nTCN node data loaded.')
# %% creating node-based annotated data object
nodeIds = list(set(df_nodes['unique ID'].values))
check = df_nodes.groupby('unique ID')['pullback_set_label'].nunique()
bad_ids = check[check > 1]

if len(bad_ids) > 0:
    print('\nFound multiple pullback_set_label for unique ID values. Please check.')
    sys.exit(0)

metaCols = ['Relative FoC change',
            'Only last FoC measurment for each sample',
            'Condition']

dictNodeData = {}
arr = None
PSL = []
print('\nCalculating node-based correlations...')
for node in nodeIds:
    selectNodes = df_nodes[df_nodes['unique ID'] == node]
    nodeSamples = selectNodes['node_elements'].tolist()
    adataSamples = adata[adata.obs.index.isin(nodeSamples)]
    '''node metadata, provided no bad_ids'''
    PSL.append(selectNodes['pullback_set_label'].unique()[0])
    # process only numerical data
    if len(nodeSamples) == 1:

        for meta_col in metaCols:
            if ' FoC ' in meta_col:
                # numerical
                mean_val = adataSamples.obs[meta_col].values.astype(float)[0]

            else:
                # categorical
                mean_val = '%s: 1.00' % adataSamples.obs[meta_col].values[0]

            if meta_col not in dictNodeData:
                dictNodeData[meta_col] = [mean_val]

            else:
                dictNodeData[meta_col].append(mean_val)

        ftrs_mean = adataSamples.X

    else:
        for meta_col in metaCols:
            if ' FoC ' in meta_col:
                # numerical
                mean_val = np.mean(
                    adataSamples.obs[meta_col].values.astype(float))
            else:
                # categorical
                categories = adataSamples.obs[meta_col].values.tolist()
                set_vls = list(set(categories))
                counts = [categories.count(x) for x in set_vls]
                frac_vls = list(map(lambda x, y: '%s: %.2f' % (
                    x, round(float(y)/sum(counts), 2)), set_vls, counts))
                mean_val = '; '.join(frac_vls)

            if meta_col not in dictNodeData:
                dictNodeData[meta_col] = [mean_val]

            else:
                dictNodeData[meta_col].append(mean_val)

        tmp = pd.DataFrame(adataSamples.X, columns=adata.var_names)
        tmp_mean = tmp.mean(axis=0)
        ftrs_mean = tmp_mean.values

    arr = (np.vstack((arr, ftrs_mean)) if (arr is not None) else ftrs_mean)
# %% construct anndata object
node_adata = ad.AnnData(arr)
node_adata.obs.index = nodeIds
node_adata.var_names = adata.var_names
node_adata.varm['protein_description'] = adata.varm['protein_description']
'''adding node-related metadata'''
node_adata.obs['unique ID'] = nodeIds
node_adata.obs['pullback_set_label'] = PSL
node_adata.uns['TDA network'] = tdaFl

for meta_col in metaCols:
    node_adata.obs[meta_col] = dictNodeData[meta_col]
# %% add transformed categorical variables
conditions = list(set(adata.obs['Condition'].values))
cat_vls = node_adata.obs['Condition'].values

for cond in conditions:
    cond_vls = list(map(lambda x: cat2num(x, cond), cat_vls))
    node_adata.obs['Enrichment with \'%s\'' % cond] = cond_vls
# %% calculate correlations and store them as varm layers
meta_layers = ['Relative FoC change', 'Only last FoC measurment for each sample', "Enrichment with 'Epinephrine 1uM'",
               "Enrichment with 'DMSO'", "Enrichment with 'Not treated'", "Enrichment with 'Doxorubicin 5 uM'"]

for meta in meta_layers:
    meta_vls = node_adata.obs[meta].values
    df_nodes = pd.DataFrame(node_adata.X, columns=node_adata.var_names)
    '''Pearson'''
    corrPearson = df_nodes.apply(lambda x: corr_fun(
        calc_type='Pearson', meta_vls=meta_vls, vls=x), axis=0)
    node_adata.varm['Pearson: %s' % meta] = corrPearson.iloc[0].values
    node_adata.varm['Pearson p-value: %s' % meta] = corrPearson.iloc[1].values
    '''Spearman'''
    corrSpearman = df_nodes.apply(lambda x: corr_fun(
        calc_type='Spearman', meta_vls=meta_vls, vls=x), axis=0)
    node_adata.varm['Spearman: %s' % meta] = corrSpearman.iloc[0].values
    node_adata.varm['Spearman p-value: %s' %
                    meta] = corrSpearman.iloc[1].values

timeExec = round((timer() - startTime), 3)
print('\n Node-based analysis completed in %.2f seconds' % timeExec)
# %% save processed data
saveOutput = True

if saveOutput:
    fl = os.path.join(rootPath, 'AnnData_%s.h5ad' % (nodeFl.split('.csv')[0]))
    node_adata.write(filename=fl)
    print('\nAnnData are saved in %s' % fl)
    '''save node data as csv'''
    node_df = node_adata.to_df()          # proteins as columns, EHTs as rows
    node_df = node_adata.obs.join(node_df)     # append obs values
    fl = os.path.join(rootPath, "Supp_TNE.csv")
    node_df.to_csv(fl, index=False)
    '''save node-based correlations as csv'''
    varm_df = pd.DataFrame({"Feature": node_adata.var_names})
    # Add all varm arrays as columns
    for key, value in node_adata.varm.items():
        varm_df[key] = value.ravel()   # or value.squeeze()

    fl = os.path.join(rootPath, "Supp_TNE_ProteinCorrelations.csv")
    varm_df.to_csv(fl, index=False)

# %% output csv: for each condition keep only features with p<0.05
saveCSV = False

if saveCSV:
    p_thresh = 5e-2
    for meta in meta_layers:
        print('\nPreparing correlations for %s...' % meta)
        '''prepare output template'''
        for corrType in ['Spearman', 'Pearson']:
            meta_output = pd.DataFrame(
                adata.var_names.tolist(), columns=['Feature'])
            for varm_val in ['gene_id', 'protein_accession', 'protein_description', 'protein_group_id', 'protein_name']:
                meta_output[varm_val] = adata.varm[varm_val]

            meta_output['%s %s' % (
                corrType, meta)] = node_adata.varm['%s: %s' % (corrType, meta)]
            meta_output['%s p-value: %s' %
                        (corrType, meta)] = node_adata.varm['%s p-value: %s' % (corrType, meta)]
            '''Step 1: removing features with p > threshold'''
            meta_output_pthresh = meta_output[meta_output['%s p-value: %s' % (
                corrType, meta)] < p_thresh].copy()

            print('\nApplied filter of p-value %s.' % str(p_thresh))
            print('\n%i features pass the threshold.' %
                  len(meta_output_pthresh))
            # add an extra column for Reactome
            meta_output_pthresh['first_protein_accession'] = (
                meta_output_pthresh['Feature'].apply(get_single_protein_acc))
            # save the output
            flOut = os.path.join(
                rootPath, '%s_NodeBasedCorr_%s.csv' % (corrType, meta))
            meta_output_pthresh.to_csv(flOut, index=False)

            print('\nOutput saved in %s.' % flOut)
# %% compare FoC correlations
if performRawCorr:
    var1 = node_adata.varm['Spearman: Only last FoC measurment for each sample']
    var2 = adata.varm['Spearman: Only last FoC measurment for each sample']

    pears_corr = np.corrcoef(var1, var2)
    print(pears_corr)

    spearm_corr = stats.spearmanr(var1, var2)
    print('Correlation: %f' % spearm_corr.correlation)
    print('p-value: %f' % spearm_corr.pvalue)

print('\nDone.')

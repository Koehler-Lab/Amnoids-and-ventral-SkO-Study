# Extended Data Fig. 10 | The role of LPM in amnioid recovery during mechanical perturbation and in cvSkOs’ mechanosensitive gene signatures.
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp
import scanpy as sc
gp.__version__


labels=["Skin-EpiVsAE"]
labels1=["AEgroup"]       
labels2=["SkinEpiderm"]

j=0
label=labels[j]
label1=labels1[j]
label2=labels2[j]
label
label1
label2

# Extended Figure 10e. Gene set enrichment analysis (GSEA) of ion channel expression in amnion (AE, AE/SE1/2) versus epidermal clusters (earlyKC, Basal KC1/2, Periderm) across integrated single-cell datasets (day 6, 28, 100+). 
# Over-representation analysis by enrichR.
# Amnion clusters show elevated expression of Piezo1/2 and calcium signaling channels (Orai, Ryanodine receptors).  
adata = sc.read_h5ad(label + ".h5ad") # data from SeuratData::ifnb
adata.obs.head()
adata.X 
adata.layers['counts'] = adata.X # Save raw counts
adata.layers['lognorm'] = adata.X
adata.obs.groupby('orig.ident')['type'].value_counts()
adata.obs.groupby('orig.ident')['CTepi'].value_counts()
adata.obs['stim'] = pd.Categorical(adata.obs['CTepi'], categories=[label2,label1], ordered=True)
adata.obs.groupby('orig.ident')['stim'].value_counts()
indices = adata.obs.sort_values(['stim']).index
adata = adata[indices,:]
bdata=adata.copy()


ions = gp.read_gmt(path="ionchannels_fam.gmt")
gene_sets=[ions]
setslabel=["IonChannel"]

jj=0
set=gene_sets[jj]
setlabel=setslabel[jj]
set
setlabel

# DEG Analysis
sc.tl.rank_genes_groups(bdata,
                        groupby='stim',
                        use_raw=False,
                        layer='lognorm',
                        method='wilcoxon',
                        groups=[label2],
                        reference=label1)

bdata.X.max() # already log1p
sc.pl.rank_genes_groups(bdata, n_genes=25, sharey=False,save=setlabel+'_'+label+'.png')

# get deg result
result = bdata.uns['rank_genes_groups']
groups = result['names'].dtype.names
degs = pd.DataFrame(
    {group + '_' + key: result[key][group]
    for group in groups for key in ['names','scores', 'pvals','pvals_adj','logfoldchanges']})

degs.head()
degs.shape

# Over-representation analysis using Enrichr
# subset up or down regulated genes
degs_sig = degs[degs[degs.columns[3]] < 0.05]
degs_up = degs_sig[degs_sig[degs_sig.columns[4]] > 0]
degs_dw = degs_sig[degs_sig[degs_sig.columns[4]] < 0]

degs_up.shape
degs_dw.shape

# Enricr API
enr_up = gp.enrichr(degs_up[degs_up.columns[0]],
                    gene_sets=set, # background=bdata.var_names
                    cutoff=0.5,
                    outdir=setlabel+'/'+label+'/Up')
enr_up.res2d.Term = enr_up.res2d.Term.str.split(" \(GO").str[0]


enr_dw = gp.enrichr(degs_dw[degs_dw.columns[0]],
                    gene_sets=set,
                    cutoff=0.5,
                    outdir=setlabel+'/'+label+'/Down')
enr_dw.res2d.Term = enr_dw.res2d.Term.str.split(" \(GO").str[0]


# concat results
enr_up.res2d['UP_DW'] = "UP"
enr_dw.res2d['UP_DW'] = "DOWN"
enr_up.res2d[enr_up.res2d.columns[3:9]]
enr_up.res2d.sort_values(by=['Adjusted P-value'])[enr_up.res2d.columns[3:9]]
enr_res = pd.concat([enr_up.res2d.sort_values(by=['Adjusted P-value']).head(), enr_dw.res2d.sort_values(by=['Adjusted P-value']).head()])
enr_up.res2d.shape
enr_dw.res2d.shape
enr_res.shape
enr_up.res2d.to_csv(setlabel+'/'+label+'/Up_0.5.csv', index=False)
enr_dw.res2d.to_csv(setlabel+'/'+label+'/Down_0.5.csv', index=False)
enr_res.to_csv(setlabel+'/'+label+'/'+label+setlabel+'_UpDown_0.5_top.csv', index=False)

from gseapy import dotplot
from gseapy.scipalette import SciPalette
sci = SciPalette()
NbDr = sci.create_colormap()

ax = gp.dotplot(enr_res,figsize=(3,5),
                x='UP_DW',
                x_order = ["UP","DOWN"],
                title=setlabel,
                cmap = NbDr.reversed(),
                cutoff=0.5,
                show_ring=True,ofname=setlabel+'/'+label+'/dotplot_0.5.png')

ax = gp.dotplot(enr_res,figsize=(3,5),
                x='UP_DW',
                x_order = ["UP","DOWN"],
                title=setlabel,
                cmap = NbDr.reversed(),
                cutoff=0.5,
                size=3,
                ofname=setlabel+'/'+label+'/dotplot2_0.5.png')


import session_info
session_info.show()
-----
gseapy              1.1.3
matplotlib          3.8.2
networkx            3.2.1
numpy               1.26.3
pandas              2.1.4
scanpy              1.10.1
session_info        1.0.0
-----
Python 3.11.5 | packaged by conda-forge | (main, Aug 27 2023, 03:34:09) [GCC 12.3.0]
Linux-4.18.0-553.63.1.el8_10.x86_64-x86_64-with-glibc2.28
-----


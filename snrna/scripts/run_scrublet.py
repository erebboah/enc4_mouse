import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

input_dir = '../scrublet/'

batches = np.arange(41,57,1)
batches = batches.astype(str)
 
for batch in batches:
    counts_matrix = scipy.io.mmread(input_dir + batch + "_matrix.mtx").T.tocsc()
    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, 
                                                                  min_cells=1, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=30)
    print(batch)           
    meta = pd.read_csv(input_dir + batch + "_barcodes.tsv",sep='\t',header = None)
    meta['doublet_scores'] = doublet_scores
    meta['doublets'] = predicted_doublets
    meta.to_csv(input_dir + batch + "_barcodes_scrublet.tsv", sep='\t', index = False,header = None)

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

meta = pd.read_csv("../ref/enc4_mouse_snrna_metadata.tsv",sep='\t')

input_dir = '../scrublet/'
for filename in os.listdir(input_dir):
    if filename.endswith((".mtx")):
        counts_matrix = scipy.io.mmread(input_dir + filename).T.tocsc()
        scrub = scr.Scrublet(counts_matrix)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, min_cells=1, min_gene_variability_pctl=85, n_prin_comps=30)
        meta = pd.read_csv(input_dir + filename[:-4] + "_metadata.csv")
        meta['doublet_scores'] = doublet_scores
        meta['doublets'] = predicted_doublets
        meta.to_csv(input_dir + filename[:-4] + "_metadata.csv", sep=',', index = False)

# Mouse snRNA Integrative Analysis
## Adrenal
### Data
- [Adrenal data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/adrenal_minimal_metadata.tsv)

### Aims
[integrate_parse_10x.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_parse_10x.R):
1. Read in pre-processed Parse and 10x data and merge counts matrices across experiments (within the same technology) for each tissue.
2. Filter nuclei by # genes, # UMIs, percent mitochondrial gene expression, and doublet score. See [detailed metadata](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs. The adrenal gland filters are the same as the brain filters. **Also filter 10x nuclei for those passing snATAC filters.**
3. Run SCT on the 3 objects to regress `percent.mt` and `nFeature_RNA`. Use  `method = "glmGamPoi"` to speed up this step, and save pre-integrated data in `seurat` folder.
4. Combine Parse standard, Parse deep, and 10x data by CCA integration. Use Parse standard as reference dataset because it contains all timepoints, while 10x data only contains 2 timepoints. 
5. Score nuclei by cell cycle using these [mouse cell cycle genes](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/mouse_cellcycle_genes.rda) to aid in manual celltype annotation.

In [adrenal gland notebook](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/ADR_snRNA.ipynb): Manual celltype annotation. 

### Results 
- `gen_celltype` and `celltypes` are basically the same, such as "Cortex" and "Medulla" (with the exception of "Myeloid" vs. "Macrophages" and "Myonuclei" vs. "Skeletal_muscle")
- `subtypes` breaks down the cortical cells into "Cortex_ZG" (zona glomerulosa) and "Cortex_ZF" (zona fasciculata) using ZG marker gene Cyp11b2. Cyp11a1 and Cyp11b1 are expressed in both layers.
- "Medulla" breaks down into "Medulla_NE" (norepinephrine-producing) and "Medulla_EPI" (epinephrine-producing) using expression of Pnmt, the enzyme responsible for concerting norepinephrine to epinephrine, which is expressed only in the epinephrine-producing cells.
- Cross-referenced [cell type marker genes](https://panglaodb.se/markers.html) with [cluster markers](https://github.com/erebboah/enc4_mouse/blob/master/snrna/seurat/markers/adrenal_cluster_marker_genes_only.pos_min.pct0.25_logfc.threshold0.25.tsv) to figure out clusters such as hepatocytes and fibroblasts.
- Cell cycle scoring helps assign cycling clusters.
- We noticed predominantly female and male cortex that we labeled "X_zone" and "Male_only_ZF", respectively, with the X-zone already described in literature.

## Left cortex
### Data
- [Cortex data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/cortex_minimal_metadata.tsv)

### Aims
[integrate_parse_10x.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_parse_10x.R):
1. Read in pre-processed Parse and 10x data and merge counts matrices across experiments (within the same technology) for each tissue.
2. Filter nuclei by # genes, # UMIs, percent mitochondrial gene expression, and doublet score. See [detailed metadata](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs.
3. Run SCT on the 3 objects to regress `percent.mt` and `nFeature_RNA`. Use  `method = "glmGamPoi"` to speed up this step, and save pre-integrated data in `seurat` folder. **Also filter 10x nuclei for those passing snATAC filters.**
4. Combine Parse standard, Parse deep, and 10x data by CCA integration. Use Parse standard as reference dataset because it contains all timepoints, while 10x data only contains 2 timepoints. 
5. Score nuclei by cell cycle using these [mouse cell cycle genes](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/mouse_cellcycle_genes.rda) to aid in manual celltype annotation.

[predict_cortex_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_cortex_celltypes.R): Use an [external 10x brain atlas](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) to predict celltype labels. The 1.1M cell dataset was subset for 1,000 cells in each celltype for a ~250,000 cell dataset (code coming soon).

In [left cortex notebook](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/CX_snRNA.ipynb): Manual celltype annotation by assigning each cluster to the celltype predicted for the majority of cells in the cluster, then adjusting the labels as we see fit.

### Results
- We decided on 3 levels of annotation: `gen_celltype` or general celltype (e.g. "Neuron"), `celltypes` for higher resolution (e.g. "Inhibitory"), and finally `subtypes` for the highest resolution of celltype annotations (e.g. "Pvalb"). 
- The external atlas did not separate their oligodendrocytes into OPCs, MFOLs, and MOLs, but we use marker genes reported in literature and on [mousebrain.org](http://www.mousebrain.org/adolescent/celltypes.html) to check marker genes and assign cell type labels.

## Hippocampus
### Data
- [Hippocampus data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/hippocampus_minimal_metadata.tsv)

### Aims
[integrate_parse_10x.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_parse_10x.R):
1. Read in pre-processed Parse and 10x data and merge counts matrices across experiments (within the same technology) for each tissue.
2. Filter nuclei by # genes, # UMIs, percent mitochondrial gene expression, and doublet score. See [detailed metadata](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs. **Also filter 10x nuclei for those passing snATAC filters.**
3. Run SCT on the 3 objects to regress `percent.mt` and `nFeature_RNA`. Use  `method = "glmGamPoi"` to speed up this step, and save pre-integrated data in `seurat` folder.
4. Combine Parse standard, Parse deep, and 10x data by CCA integration. Use Parse standard as reference dataset because it contains all timepoints, while 10x data only contains 2 timepoints. 
5. Score nuclei by cell cycle using these [mouse cell cycle genes](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/mouse_cellcycle_genes.rda) to aid in manual celltype annotation.

[predict_hippocampus_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_hippocampus_celltypes.R): Use an [external 10x brain atlas](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) to predict celltype labels. The 1.1M cell dataset was subset for 1,000 cells in each celltype for a ~250,000 cell dataset (code coming soon).

In [hippocampus notebook](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/HC_snRNA.ipynb): Manual celltype annotation by assigning each cluster to the celltype predicted for the majority of cells in the cluster, then adjusting the labels as we see fit.

### Results
- We decided on 3 levels of annotation: `gen_celltype` or general celltype (e.g. "Neuron"), `celltypes` for higher resolution (e.g. "Inhibitory"), and finally `subtypes` for the highest resolution of celltype annotations (e.g. "Pvalb"). 
- The external atlas did not separate their oligodendrocytes into OPCs, MFOLs, and MOLs, but we use marker genes reported in literature and on [mousebrain.org](http://www.mousebrain.org/adolescent/celltypes.html) to check marker genes and assign cell type labels.

## Heart
### Data
- [Heart data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/cortex_minimal_metadata.tsv)

### Aims
[integrate_parse_10x.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_parse_10x.R):
1. Read in pre-processed Parse and 10x data and merge counts matrices across experiments (within the same technology) for each tissue.
2. Filter nuclei by # genes, # UMIs, percent mitochondrial gene expression, and doublet score. See [detailed metadata](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs. **Also filter 10x nuclei for those passing snATAC filters.**
3. Run SCT on the 3 objects to regress `percent.mt` and `nFeature_RNA`. Use  `method = "glmGamPoi"` to speed up this step, and save pre-integrated data in `seurat` folder.
4. Combine Parse standard, Parse deep, and 10x data by CCA integration. Use Parse standard as reference dataset because it contains all timepoints, while 10x data only contains 2 timepoints. 
5. Score nuclei by cell cycle using these [mouse cell cycle genes](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/mouse_cellcycle_genes.rda) to aid in manual celltype annotation.

[predict_heart_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_cortex_celltypes.R): use 2 external datasets to predict celltypes.

- [Stressed mouse ventricles](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8810/files/) and associated [paper](https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.119.045115)
- [Human heart cell atlas](https://www.heartcellatlas.org/) converted to mouse gene IDs.

In [heart notebook](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/HT_snRNA.ipynb): Manual celltype annotation by assigning each cluster to the celltype predicted for the majority of cells in the cluster, then adjusting the labels as we see fit.


### Results
- Stressed heart data did not include annotations, so I annotated them manually in [predict_heart_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_cortex_celltypes.R) using the marker genes in [Fig. 1C](https://www.ahajournals.org/cms/asset/5d52d0eb-9c13-4a10-81f4-58e18559d64d/circulationaha.119.045115.fig01.jpg)
- Main heart celltypes are CM (cardiomyocytes), CF (cardiac fibroblasts), and CE (cardiac endothelial)
- Split CM, CF, CE clusters by timepoint (early, middle, late) when appropriate
- Split other clusters by cycling
- Stressed heart data predicted most celltypes, but I used the adipocyte and atrial cardiomyocyte labelling information from the heart atlas data

## Gastrocnemius
### Data
- [Gastroc data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/gastrocnemius_minimal_metadata.tsv)

### Aims
[integrate_parse_10x.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_parse_10x.R):
1. Read in pre-processed Parse and 10x data and merge counts matrices across experiments (within the same technology) for each tissue.
2. Filter nuclei by # genes, # UMIs, percent mitochondrial gene expression, and doublet score. See [detailed metadata](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs. **Also filter 10x nuclei for those passing snATAC filters.**
3. Run SCT on the 3 objects to regress `percent.mt` and `nFeature_RNA`. Use  `method = "glmGamPoi"` to speed up this step, and save pre-integrated data in `seurat` folder.
4. Combine Parse standard, Parse deep, and 10x data by CCA integration. Use Parse standard as reference dataset because it contains all timepoints, while 10x data only contains 2 timepoints. 
5. Score nuclei by cell cycle using these [mouse cell cycle genes](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/mouse_cellcycle_genes.rda) to aid in manual celltype annotation.

[predict_gastroc_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_gastroc_celltypes.R): Use an [external TA dataset](https://www.synapse.org/#!Synapse:syn21676145/files/) from [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7733460/pdf/41467_2020_Article_20063.pdf) to predict celltype labels.

In [gastrocnemius notebook](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/GC_snRNA.ipynb): Manual celltype annotation by assigning each cluster to the celltype predicted for the majority of cells in the cluster, then adjusting the labels as we see fit. 

### Results
- Seurat CCA integration still works ok, but this tissue has the most difference between Parse and 10x, most likely due to differences in nuclei prep that are most apparent in muscle.
- Low cluster definition, especially in myonuclei

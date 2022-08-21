# Mouse snRNA Integrative Analysis

## Aims
[integrate_parse_10x.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/integrate_parse_10x.R):
1. Read in pre-processed Parse and 10x data and merge counts matrices across experiments (within the same technology) for each tissue.
2. Filter nuclei by # genes, # UMIs, percent mitochondrial gene expression, and doublet score. See [detailed metadata](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs.
3. Run SCT on the 3 objects to regress `percent.mt` and `nFeature_RNA`. Use  `method = "glmGamPoi"` to speed up this step, and save pre-integrated data in `seurat` folder.
4. Combine Parse standard, Parse deep, and 10x data by CCA integration. Use Parse standard as reference dataset because it contains all timepoints, while 10x data only contains 2 timepoints. 
5. Score nuclei by cell cycle using these [mouse cell cycle genes](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/mouse_cellcycle_genes.rda) to aid in manual celltype annotation.


## Brain
### Data
- [Hippocampus data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/hippocampus_minimal_metadata.tsv)
- [Cortex data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/cortex_minimal_metadata.tsv)

[predict_hippocampus_celltypes.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/predict_hippocampus_celltypes.R), [predict_cortex_celltypes.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/predict_cortex_celltypes.R): Use an [external 10x brain atlas](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) to predict celltype labels. The 1.1M cell dataset was subset for 1,000 cells in each celltype for a ~250,000 cell dataset (code coming soon).

[HC_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/HC_snRNA.ipynb), [CX_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/CX_snRNA.ipynb): Manual celltype annotation by assigning each cluster to the celltype predicted for the majority of cells in the cluster, then adjusting the labels as we see fit. Find marker genes for gen_celltype, celltypes, and subtypes and save in seurat/markers


### Results
- We decided on 3 levels of annotation: `gen_celltype` or general celltype (e.g. "Neuron"), `celltypes` for higher resolution (e.g. "Inhibitory"), and finally `subtypes` for the highest resolution of celltype annotations (e.g. "Pvalb"). 
- The external atlas did not separate their oligodendrocytes into OPCs, MFOLs, and MOLs, but we use marker genes reported in literature and on [mousebrain.org](http://www.mousebrain.org/adolescent/celltypes.html) to check marker genes and assign cell type labels.

## Adrenal
### Data
- [Adrenal data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/adrenal_minimal_metadata.tsv)

[ADR_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/ADR_snRNA.ipynb):
Manual celltype annotation. Find marker genes for `gen_celltype`, `celltypes`, and `subtypes` and save in `seurat/markers`.

### Results 
- `gen_celltype` and `celltypes` are basically the same, such as "Cortex" and "Medulla" (with the exception of "Myeloid" vs. "Macrophages" and "Myonuclei" vs. "Skeletal_muscle")
- `subtypes` breaks down the cortical cells into groups such as "Cortex_ZG" (zona glomerulosa), "Cortex_ZF" (zona fasciculata), "X_zone", and "Y_zone". "Medulla" breaks down into "Medulla_NE" (norepinephrine-producing) and "Medulla_EPI" (epinephrine-producing).
- Cross-referenced [cell type marker genes](https://panglaodb.se/markers.html) with [cluster markers](https://github.com/erebboah/enc4_mouse/blob/master/snrna/seurat/markers/adrenal_cluster_marker_genes_only.pos_min.pct0.25_logfc.threshold0.25.tsv) to figure out clusters such as hepatocytes and fibroblasts.
- Cell cycle scoring helps assign cycling clusters.
- We noticed predominantly male and female cortex that we labeled Y_zone and X_zone, respectively, with the X-zone already described in literature



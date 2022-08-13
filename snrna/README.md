# Mouse snRNA Integrative Analysis
## Hippocampus
### Data
- [Hippocampus data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/hippocampus_minimal_metadata.tsv)

### Aims
[integrate_hippocampus.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/integrate_hippocampus.R):
1. Read in pre-processed Parse and 10x data and merge counts matrices across experiments (within the same technology) for each tissue.
2. Filter nuclei by # genes, # UMIs, percent mitochondrial gene expression, and doublet score. See [detailed metadata](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs.
3. Run SCT on the 3 objects to regress `percent.mt` and `nFeature_RNA`. Use  `method = "glmGamPoi"` to speed up this step, and save pre-integrated data in `seurat` folder.
4. Combine Parse standard, Parse deep, and 10x data by CCA integration. Use Parse standard Seurat object as reference dataset because it contains all timepoints, while 10x data only contains 2 timepoints. 

[predict_hippocampus_celltypes.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/predict_hippocampus_celltypes.R)
1. Use an [external 10x brain atlas](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) to predict celltype labels. The 1.1M cell dataset was subsetted by 1,000 cells in each celltype for a ~250,000 cell dataset (code coming soon).
2. Score nuclei by cell cycle using these [mouse cell cycle genes](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/mouse_cellcycle_genes.rda) to aid in manual celltype annotation.

**In this notebook**:
Manual celltype annotation by assigning each cluster to the celltype predicted for the majority of cells in the cluster, then adjusting the labels as we see fit.


### Results
- Seurat CCA works pretty well for integrating the 3 types of experiments: Parse standard, Parse deep, and 10x multiome. 
- We decided on 3 levels of annotation: `gen_celltypes` or general celltypes (e.g. "Neuron"), `celltypes` for higher resolution (e.g. "Inhibitory"), and finally `subtypes` for the highest resolution of celltype annotations (e.g. "Pvalb"). 
- The external atlas did not separate their oligodendrocytes into OPCs, MFOLs, and MOLs, but we use our expertise with the brain to check marker genes and assign cell type labels.

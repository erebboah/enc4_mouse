# Mouse snRNA Integrative Analysis
## Hippocampus
### Data
- [Hippocampus data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/hippocampus_minimal_metadata.tsv)

### Aims
1. Reads in pre-processed Parse and 10x data and merge counts matrices across experiments (within the same technology) for each tissue.
2. Combine Parse standard, Parse deep, and 10x data by CCA integration.
3. Use an [external 10x brain atlas](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) to predict celltype labels.
4. Manual celltype annotation by assigning each cluster to the celltype predicted for the majority of cells in the cluster, then adjusting the labels as we see fit.

### Results
- Seurat CCA works pretty well for integrating the 3 types of experiments: Parse standard, Parse deep, and 10x multiome. 
- We decided on 3 levels of annotation: `gen_celltypes` or general celltypes (e.g. "Neuron"), `celltypes` for higher resolution (e.g. "Inhibitory"), and finally `subtypes` for the highest resolution of celltype annotations (e.g. "Pvalb"). 
- The external atlas did not separate their oligodendrocytes into OPCs, MFOLs, and MOLs, but we use our expertise with the brain to check marker genes and assign cell type labels.

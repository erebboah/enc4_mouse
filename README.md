# ENCODE 4 Mouse snRNA-seq and snATAC-seq Analysis
## Overview
Analysis of postnatal timecourse of Bl6/Cast F1 hybrid mouse development in 5 core tissues: cortex, hippocampus, heart, adrenal gland, and gastrocnemius, at 7 timepoints for snRNA-seq: PND04, PND10, PND14, PND25, PND36, PNM02, and PNM18-20 and 2 timepoints for multiome snRNA-seq+snATAC-seq: PND14 and PNM02, with 2 males and 2 females per timepoint. Also added C2C12 mouse myoblast cell line data at 0hr and 72hr differentiation timepoints.

## Contact info
Analyst: Liz Rebboah, Mortazavi Lab. erebboah@uci.edu

## snRNA-seq
### Data
- ENCODE carts: [Parse](https://www.encodeproject.org/carts/enc4_mouse_snrna_parse/), [10x](https://www.encodeproject.org/carts/enc4_mouse_snrna_10x/)
- [Hippocampus data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/hippocampus_minimal_metadata.tsv)
- [Cortex data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/cortex_minimal_metadata.tsv)
- [Adrenal data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/adrenal_minimal_metadata.tsv)
- [Heart data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/heart_minimal_metadata.tsv)
- [Gastrocnemius data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/gastrocnemius_minimal_metadata.tsv)
- [C2C12 data table](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/c2c12_minimal_metadata.tsv)

Note: All Parse snRNA-seq experiments were done at the Mortazavi lab while all 10x multiome experiments were done at the Snyder lab, from nuclei isolation to sequencing.

### Pre-processing
1. [step1_get_data.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/step1_get_data.sh) fetches **unfiltered sparse gene count matrix of all reads** tar.gz files from the ENCODE portal using the batch download script `xargs -L 1 curl -O -J -L < files.txt` for [Parse](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/encode_parse_scrna_files.txt) and [10x](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/encode_10x_scrna_files.txt) separately. The tar folders are unzipped in `counts_parse` or `counts_10x` with the ENCODE "ENCFF" filename ID as the folder name. There are 40 10x folders and 436 Parse folders. **Unzipping the data requires 39G of space for Parse, 32G for 10x.** Next [get_counts_parse.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/get_counts_parse.R) is called to merge Parse data by experimental batch across file IDs and filter the resulting counts matrices by gene biotypes and for nuclei > 500 UMI, which are saved as sparse matrices and tsv files for genes and barcodes in the `scrublet` folder.  Finally [get_counts_10x.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/get_counts_10x.R) reads in 10x data by file ID and filters the counts matrices by gene biotypes and for nuclei > 0 UMI for ambient RNA removal. The sparse matrices and tsv files for genes and barcodes are saved in the same file ID folder in `counts_10x`.  
2. [step2_cellbender.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/step2_cellbender.sh) runs Cellbender as an array job to remove ambient RNA from the 40 10x datasets using [these Cellbender settings](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/10x_cellbender_settings.csv), where the first column is path to the file ID folder, second column is expected cells, and third column is total droplets included. The unfiltered Cellbender output is saved as cellbender.h5 in the same file ID folder in `counts_10x`. **The Cellbender output requires an additional 56G of space.**
3. [step3_run_scrublet.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/step3_run_scrublet.sh) runs followup script [format_cellbender_output.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/format_cellbender_output.R) to output >500 UMI sparse matrices to the `scrublet` folder from the 10x cellbender.h5 files, matching the status of the Parse data. The python script [run_scrublet.py](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/run_scrublet.py) detects doublets (same script for Parse and 10x) and outputs a modified barcodes tsv (`_barcodes_scrublet.tsv`) file in the `scrublet` folder with additonal columns for doublet information.
4. [step4_submit_data.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/step4_submit_data.sh) runs [Submission_files.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/Submission_files.R) to output files for submission to Synapse/ENCODE portal using [these guidelines](https://docs.google.com/document/d/1NAuihlO3Qd4y588Fpv1iWTyO2jxo2u8M9U-vfwnGbmE/edit#).


### Celltype Annotation
**Note: must run snATAC script [step1_atac_archr.sh](https://github.com/erebboah/enc4_mouse/blob/master/snatac/scripts/step1_get_data.sh) before intergrating data so we can filter the 10x snRNA nuclei for those also passing snATAC filters.**

Each tissue is integrated across technologies and annotated using Seurat. The final annotations have 3 levels of granularity: `gen_celltype`, `celltypes`, and `subtypes`. See [celltype metadata](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/enc4_mouse_snrna_celltypes_c2c12.csv) to see how these levels relate to each other.

[integrate_parse_10x.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_parse_10x.R) merges counts across technologies and makes 3 Seurat objects for Parse standard, Parse deep, and 10x. Nuclei are filtered (see [detailed metadata](https://github.com/erebboah/enc4_mouse/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs). CCA integrates the 3 objects. The `combined.sct` object is processed with PCA, UMAP, SNN graph construction, and high-resolution clustering. Marker genes are called from Seurat clusters and saved in the `seurat` folder.

#### [Hippocampus](https://www.synapse.org/#!Synapse:syn39422435)
1. Run integration R script with [integrate_hippocampus.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_hippocampus.sh).
2. Check integration results and clustering resolution in [HC_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/HC_snRNA.ipynb).
3. [predict_hippocampus_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_hippocampus_celltypes.R) uses an [external 10x dataset](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) from mouse hippcampus and cortex [subsetted](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/subsample_brain_data.R) by 1,000 nuclei in each annotated subtype to predict celltypes.
4. Check prediction results in [HC_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/HC_snRNA.ipynb) and make adjustments.

#### [Cortex](https://www.synapse.org/#!Synapse:syn39422426)
1. Run integration R script with [integrate_cortex.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_cortex.sh). 
2. Check integration results and clustering resolution in [CX_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/CX_snRNA.ipynb).
3. [predict_cortex_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_cortex_celltypes.R) uses an [external 10x dataset](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) from mouse hippcampus and cortex [subsetted](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/subsample_brain_data.R) by 1,000 nuclei in each annotated subtype to predict celltypes.
4. Check prediction results in [CX_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/CX_snRNA.ipynb) and make adjustments.

#### [Adrenal](https://www.synapse.org/#!Synapse:syn39422421)
1. Run integration R script with [integrate_adrenal.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_adrenal.sh). 
2. Check integration results and clustering resolution in [ADR_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/ADR_snRNA.ipynb).
3. Manually annotate `gen_celltype`, `celltypes`, and `subtypes` in [ADR_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/ADR_snRNA.ipynb).

#### [Heart](https://www.synapse.org/#!Synapse:syn39422439)
1. Run integration R script with [integrate_heart.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_heart.sh). 
2. Check integration results and clustering resolution in [HT_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/HT_snRNA.ipynb).
3. [predict_heart_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_heart_celltypes.R) uses 2 external datasets to predict celltypes:
    - [Stressed mouse ventricles](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8810/files/) and associated [paper](https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.119.045115)
    - [Human heart cell atlas](https://www.heartcellatlas.org/) converted to mouse.
5. Check prediction results in [HT_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/HT_snRNA.ipynb) and make adjustments.

#### [Gastrocnemius](https://www.synapse.org/#!Synapse:syn39422445)
1. Run integration R script with [integrate_gastroc.sh](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/integrate_gastroc.sh). 
2. Check integration results and clustering resolution in [GC_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/GC_snRNA.ipynb).
3. [predict_gastroc_celltypes.R](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/predict_gastroc_celltypes.R) uses external dataset from TA muscle to predict celltypes:
    - [PND10, PND21, PNM 05](https://www.synapse.org/#!Synapse:syn21676145/files/) and associated [paper](https://www.nature.com/articles/s41467-020-20063-w5)
5. Check prediction results in [GC_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/GC_snRNA.ipynb) and make adjustments.

#### C2C12
1. Add metadata from paper in [C2C12_snRNA.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/C2C12_snRNA.ipynb).

### Regulatory gene list
1. Correlate miRNA expression to host gene expression: [miRNA_hostGene_correlation.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/miRNA_hostGene_correlation.ipynb)
2. Download and concatenate GO terms: [Histone_regulators.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/Histone_regulators.ipynb)
3. Make final regulatory gene "display table": [Regulatory_genes.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/Regulatory_genes.ipynb)
4. Subset data on regulatory genes and merge across tissues: [snRNA_regulatory_subset.ipynb](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/snRNA_regulatory_subset.ipynb)

## snATAC-seq
1. [step1_atac_archr.sh](https://github.com/erebboah/enc4_mouse/blob/master/snatac/scripts/step1_get_data.sh) downloads the 10x multiome fragment [files](https://github.com/erebboah/enc4_mouse/blob/master/snatac/ref/encode_10x_fragment_files.txt) from the ENCODE portal using the batch download script `xargs -L 1 curl -O -J -L < files.txt` from this [cart](https://www.encodeproject.org/carts/enc4_mouse_snatac/) for the 5 tissues. The tar folders are unzipped in the `fragments` folder with the ENCODE "ENCFF" filename ID as the folder name. There are 40 folders. **Unzipping the data requires 93G of space.**. Calls [make_archr_proj.R](https://github.com/erebboah/enc4_mouse/blob/master/snatac/scripts/make_archr_proj.R) which makes ArchR Arrow files, initializes ArchR Project and filters nuclei (minTSS = 4, minFrags = 1000, filterDoublets). 
2. [step2_archr_process.sh](https://github.com/erebboah/enc4_mouse/blob/master/snatac/scripts/step2_archr_process.sh) calls [process_archr_proj.R](https://github.com/erebboah/enc4_mouse/blob/master/snatac/scripts/process_archr_proj.R) to also filter snATAC nuclei for those that are also in the corresponding [10x snRNA metadata](https://github.com/erebboah/enc4_mouse/blob/master/snrna/seurat/all_tissues_10x_TFs_mirhgs_chromreg_metadata.csv) (passed RNA filters). Also includes standard ArchR processing such as iterative LSI, clustering, and UMAPs. **All snATAC data is multiome, so cell type annotations are simply transferred from the snRNA annotation to the ArchR object.**
3. [step3_tissue_analysis.sh](https://github.com/erebboah/enc4_mouse/blob/master/snatac/scripts/step3_tissue_analysis.sh) calls [archr_tissue_objects.R](https://github.com/erebboah/enc4_mouse/blob/master/snatac/scripts/archr_tissue_objects.R) to split 10x mouse ArchR object into 5 separate tissue objects followed by standard ArchR processing.

## Submission files
[Code for submission files](https://github.com/erebboah/enc4_mouse/blob/master/snrna/scripts/Submission_files.R). Files are currently hosted on Synapse:
- [Hippocampus](https://www.synapse.org/#!Synapse:syn39422435)
- [Cortex](https://www.synapse.org/#!Synapse:syn39422426)
- [Adrenal](https://www.synapse.org/#!Synapse:syn39422421)
- [Heart](https://www.synapse.org/#!Synapse:syn39422439)
- [Gastrocnemius](https://www.synapse.org/#!Synapse:syn39422445)

## Figure captions
Text captions for figures included in figures tarball. Not all tissues have the same figures in `figures/snrna/annotation`.

### figures/snatac/
- `Plot-UMAP-Sample-Clusters.pdf`: 4 snATAC UMAPs grouped by metadata: sample ID, tissue, timepoint, and annotated general celltype, celltypes, and subtypes.
### figures/snrna/annotation
- `UMAP_final_celltypes.pdf`: snRNA UMAP grouped by final annotated celltypes.
- `UMAP_final_gen_celltype.pdf`: snRNA UMAP grouped by final annotated general celltype.
- `UMAP_final_subtypes.pdf`: snRNA UMAP grouped by final annotated general subtypes.
- `UMAP_maximum_predictions.pdf`: snRNA UMAP grouped by the maximum predicted celltype in each cluster.
- `UMAP_predictions.pdf`: snRNA UMAP grouped by predicted celltypes, if using an external reference.
- `celltypes_marker_dotplot.pdf`: Dot plot of cells grouped by annotated celltype showing average expression and percent of cells expressing selected marker genes.
- `cluster_marker_dotplot.pdf`: Dot plot of cells grouped by cluster showing average expression and percent of cells expressing selected marker genes.
- `marker_featureplots.pdf`: snRNA UMAP "feature plots" of selected marker genes.
- `subtypes_marker_dotplot.pdf`: Dot plot of cells grouped by annotated subtypes showing average expression and percent of cells expressing selected marker genes.
- `timepoint_celltypes_proportions.pdf`: Proportion of celltypes across postnatal development.
### figures/snrna/clustering
- `Parse_10x_experiment_distribution.pdf`: snRNA UMAPs grouped by cluster and split by experiment (Parse standard/deep, 10x) to check integration results, also cluster-level experiment proportion barplot.
- `UMAP_Parse_10x.pdf`: snRNA UMAP grouped by technology and cluster.
- `UMAP_cluster_sample_barplot.pdf`: snRNA UMAP grouped by cluster and cluster-level sample proportion barplot.
- `age_sex_barplot.pdf`: snRNA UMAP grouped by developmental timepoint and sex alongside cluster-level timepoint and sex proportion barplots. 
### figures/snrna/qc
- `experiment_kneeplots.pdf`: Knee plots showing total # UMIs of ranked cells and red dashed line at 500 UMIs, grouped by snRNA experiment (label is counts matrix file accession). 
- `experiment_violinplots.pdf`: Violin plots of # genes, # UMIs, percent mitochondrial, and percent ribosomal gene expression, grouped by sample and split by experiment (Parse standard/deep, 10x).
- `qc_featureplot.png`: snRNA UMAP "feature plots" of # genes, # UMIs, percent mitochondrial and percent ribosomal gene expression, doublet score, and cell cycle G2M score.

## Auxiliary data - available upon request
Auxiliary data tarball includes processed Seurat snRNA objects before and after integration and Seurat object metadata as a csv file, Seurat cluster marker genes, and snATAC ArchR project folders (arrow files, rds files, and plots).
- Adrenal gland: 20GB
- Left cortex: 24GB
- Hippocampus: 19GB
- Heart: 22GB
- Gastrocnemius: 17GB

## Directory structure
<img src="https://github.com/erebboah/enc4_mouse/blob/master/directory_struct.png" width="355" height="590">

## FAQ
**Q:** Why is the snRNA analysis so much longer than the snATAC analysis?

**A:** The snRNA comes from 2 different single cell platforms and must be integrated, then celltypes must be manually annotated. The ATAC comes from the [10x multiome](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression) where chromatin accessibility and gene expression are detected simultaneously from the same nucleus. Therefore the same nuclei that were annotated in the snRNA analysis are in the snATAC data, so the celltype annotations were simply transferred to the ATAC Archr metadata. 


**Q:** Why are there differences in pre-processing between Parse and 10x snRNA?

**A:** The snRNA comes from [Parse Biosciences](https://www.parsebiosciences.com/) which uses combinatorial barcoding, and [10x Genomics](https://www.10xgenomics.com/), which uses droplet-based barcoding. Droplet-based barcoding introduces the possibility of RNA outside of the nucleus ending up in a droplet and getting barcoded along with real nuclei ("empty droplets"). Combinatorial barcoding requires each RNA molecule to be fixed inside the nuclei across every round of barcoding, so we feel that the empty droplets filter is unnecessary. Other than ambient RNA removal, the processing for Parse and 10x data is the same.


**Q:** Why are there so many Parse snRNA experiments?

**A:** 1. We conducted Parse snRNA-seq for all 7 postnatal timepoints, while 10x snRNA-seq was done for PND14 and PNM2. 2. We sequenced every short-read Parse experiment at a standard depth, but a subset of 1,000-2,000 nuclei were deeply sequenced. These deeply sequenced nuclei were also sequenced with either PacBio or ONT long read platforms (single-nucleus LR analysis code coming soon). See our [LR-Split-seq paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02505-w) for more details. 


**Q:** Why do you need to integrate 2 Parse Seurat RNA objects with 1 10x Parse RNA object?

**A:** See above, but basically the difference in sequencing depth between Parse "standard" and Parse "deep" is a batch effect best fixed with the same integration strategy for combining the Parse and 10x experiments. The raw counts matrices can be merged within technology, depth, and tissue (i.e. across timepoints and sexes) with no batch effects, but differences in Parse and 10x experiments (including differences in nuclei preparation) required a heavy hand at the integration step.

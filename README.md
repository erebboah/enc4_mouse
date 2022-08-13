# ENC4 Mouse snRNA-seq and snATAC-seq Analysis

## snRNA-seq
### Data
- ENCODE carts: [Parse](https://www.encodeproject.org/carts/enc4_mouse_snrna_parse/), [10x](https://www.encodeproject.org/carts/enc4_mouse_snrna_10x/)
- [Hippocampus data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/hippocampus_minimal_metadata.tsv)
- [Adrenal data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/adrenal_minimal_metadata.tsv)
- [Cortex data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/cortex_minimal_metadata.tsv)
- [Heart data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/heart_minimal_metadata.tsv)
- [Gastrocnemius data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/gastrocnemius_minimal_metadata.tsv)
- [C2C12 data table](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/c2c12_minimal_metadata.tsv)

Note: All Parse experiments were done at the Mortazavi lab while all 10x multiome experiments were done at the Snyder lab, from nuclei isolation to sequencing

### Pre-processing
1. [step1_get_data.sh](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/step1_get_data.sh) fetches **unfiltered sparse gene count matrix of all reads** tar.gz files from the ENCODE portal using the batch download script `xargs -L 1 curl -O -J -L < files.txt` for [Parse](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/encode_parse_scrna_files.txt) and [10x](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/encode_10x_scrna_files.txt) separately. The tar folders are unzipped in `counts_parse` or `counts_10x` with the ENCODE "ENCFF" filename ID as the folder name. There are 40 10x folders and 436 Parse folders. **Unzipping the data requires 39G of space for Parse, 32G for 10x.** Next [get_counts_parse.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/get_counts_parse.R) is called to merge Parse data by experimental batch across file IDs and filter the resulting counts matrices by gene biotypes and for nuclei > 500 UMI, which are saved as sparse matrices and tsv files for genes and barcodes in the `scrublet` folder.  Finally [get_counts_10x.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/get_counts_10x.R) reads in 10x data by file ID and filters the counts matrices by gene biotypes and for nuclei > 0 UMI for ambient RNA removal. The sparse matrices and tsv files for genes and barcodes are saved in the same file ID folder in `counts_10x`.  
2. [step2_cellbender.sh](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/step2_cellbender.sh) runs Cellbender as an array job to remove ambient RNA from the 40 10x datasets using [these Cellbender settings](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/10x_cellbender_settings.csv), where the first column is path to the file ID folder, second column is expected cells, and third column is total droplets included. The unfiltered Cellbender output is saved as cellbender.h5 in the same file ID folder in `counts_10x`. **The Cellbender output requires an additional 56G of space.**
3. [step3_run_scrublet.sh](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/step3_run_scrublet.sh) runs followup script [format_cellbender_output.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/format_cellbender_output.R) to output >500 UMI sparse matrices to the `scrublet` folder from the 10x cellbender.h5 files, matching the status of the Parse data. The python script [run_scrublet.py](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/run_scrublet.py) detects doublets and outputs a modified barcodes tsv (`_barcodes_scrublet.tsv`) file in the `scrublet` folder with additonal columns for doublet score and doublet True/False.

### Celltype Annotation
Each tissue is integrated across technologies and annotated using Seurat. The final annotations have 3 levels of granularity: `gen_celltype`, `celltypes`, and `subtypes`. See [celltype metadata](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/enc4_mouse_snrna_celltypes_c2c12.csv) to see how these levels relate to each other.

#### Hippocampus
1. [integrate_hippocampus.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/integrate_hippocampus.R) merges counts across technologies and makes 3 Seurat objects for Parse standard, Parse deep, and 10x. Nuclei are filtered (see [detailed metadata](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/ref/enc4_mouse_snrna_metadata.tsv) for filter cutoffs). CCA integrates the 3 objects. The `combined.sct` object is processed with PCA, UMAP, SNN graph construction, and high-resolution clustering.
2. Check integration results and clustering resolution in `HC_snRNA.ipynb`.
3. [predict_hippocampus_celltypes.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snrna/scripts/predict_hippocampus_celltypes.R) uses an [external 10x dataset](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) from mouse hippcampus and cortex subsetted by 1,000 nuclei in each annotated subtype (code coming soon) to predict celltypes. The resulting predicted.id is saved in `atlas_predictions` in the object metadata. 
4. Check prediction results in `HC_snRNA.ipynb` and make adjustments.

#### Adrenal: In progress
#### Cortex: In progress
#### Heart: In progress
#### Gastrocnemius: In progress
#### C2C12: In progress

## snATAC-seq
[step1_atac_archr.sh](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snatac/scripts/step1_atac_archr.sh) does the following:
1. Downloads the 10x multiome fragment [files](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snatac/ref/encode_10x_fragment_files.txt) from the ENCODE portal using the batch download script `xargs -L 1 curl -O -J -L < files.txt` from this [cart](https://www.encodeproject.org/carts/enc4_mouse_snatac/) for the 5 tissues. The tar folders are unzipped in the `fragments` folder with the ENCODE "ENCFF" filename ID as the folder name. There are 40 folders. **Unzipping the data requires 93G of space.** 
2. Calls [make_archr_proj.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snatac/scripts/make_archr_proj.R) which makes ArchR Arrow files, initializes ArchR Project, adds metadata, and filters nuclei. Nuclei must be present in filtered RNA object as well as passing ArchR filters (minTSS = 4, minFrags = 1000, filterDoublets). All Archr outputs, Arrow files, and project folders are in the `archr` folder.
3. With [make_archr_proj.R](https://github.com/erebboah/ENC4_Mouse_SingleCell/blob/master/snatac/scripts/make_archr_proj.R), continue processing ArchR project by dimensionality reduction, clustering, and plotting resulting UMAPs.

## Directory structure
<img src="[https://your-image-url.type](https://github.com/erebboah/enc4_mouse/blob/master/directory_struct.png)" width="100" height="100">

## FAQ
**Q:** Why are there differences in pre-processing between Parse and 10x snRNA?

**A:** Droplet-based barcoding introduces the possibility of RNA outside of the nucleus ending up in a droplet and getting barcoded along with real cells ("empty droplets"). Combinatorial barcoding requires each RNA molecule to be fixed inside the nuclei across every round of barcoding until lysis, so we feel that the empty droplets filter is unnecessary. Other than ambient RNA removal, the processing for Parse and 10x data is the same.


**Q:** Why are there so many Parse snRNA experiments?

**A:** We sequenced every short-read Parse experiment at a standard depth, but a subset of 1,000-2,000 nuclei were also sequenced deeply. These deeply sequenced nuclei were also sequenced with either PacBio or ONT long read platforms. See our [LR-Split-seq paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02505-w) for more details. 


**Q:** Why do you need to integrate 2 Parse Seurat RNA objects with 1 10x Parse RNA object?

**A:** See above, but basically the difference in sequencing depth between Parse "standard" and Parse "deep" is a batch effect best fixed with the same integration strategy for combining the Parse and 10x experiments. The raw counts matrices can be merged within technology, depth, and tissue (i.e. across timepoints and sexes) with no batch effects, but differences in Parse and 10x experiments (including differences in nuclei preparation) required a heavy hand at the integration step.


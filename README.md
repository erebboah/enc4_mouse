# ENC4_Mouse_SingleCell
``/share/crsp/lab/seyedam/share/enc4_mouse``

Scripts to import ENCODE processed data for snRNA-seq experiments and snATAC-seq experiments to Seurat and ArchR objects, respectively.

## snRNA-seq
1. ``scripts/step1_get_data.sh`` and ``scripts/merge_counts.R`` download the STARSolo unfiltered sparse gene count matrices of all reads from this [cart](https://www.encodeproject.org/carts/enc4_mouse_snrna), format with gene IDs and unique cell barcodes, merge by tissue, and **filter** by gene biotype and > 500 UMIs per nucleus.
2. ``scripts/step2_scrublet_cellbender.sh`` runs Scrublet on > 500 UMI counts matrices, split up appropriately by tissue/experiment. **Cellbender in progress**, need to integrate it within this pipeline.
3. ``scripts/step3_seurat_process.sh`` and ``scripts/seurat_process.R`` carries out final cell filtering and integrates Parse + 10x data for each tissue separately. 

## Mouse snATAC-seq
``scripts/step1_atac_archr.sh`` does the following:
1. Downloads the fragment files from this [cart](https://www.encodeproject.org/carts/enc4_mouse_snatac/) for all 5 tissues using bash syntax.
2. Calls ``scripts/make_archr_proj.R`` which makes ArchR Arrow files, initializes ArchR Project, adds metadata, and filters cells. Cells must be present in filtered RNA object as well as passing ArchR filters (minTSS = 4, minFrags = 1000, filterDoublets).
4. Continue ArchR processing by dimensionality reduction, clustering, and plotting resulting UMAPs.

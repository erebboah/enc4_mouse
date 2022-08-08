# ENC4_Mouse_SingleCell
``/share/crsp/lab/seyedam/share/enc4_mouse``

Scripts to import ENCODE processed data for snRNA-seq experiments and snATAC-seq experiments to Seurat and ArchR objects, respectively.

## snRNA-seq
### Previous steps for Parse Biosciences data
1. `step1_get_parse_data.sh` fetches Parse data from the ENCODE portal using [this cart](https://www.encodeproject.org/carts/enc4_mouse_snrna_parse/) and calls `get_counts_parse.R` to merge Parse data by experimental batch for doublet detection.
2. `run_scrublet.sh` runs python script `run_scrublet.py` to detect doublets.

## Previous steps for 10x data
1. `step1_get_10x_data.sh` fetches 10x data from the ENCODE portal using [this cart](https://www.encodeproject.org/carts/enc4_mouse_snrna_10x/) and calls `get_counts_10x.R` to output cellbender options.
2. `step2_cellbender_10x.sh` runs Cellbender to remove ambient RNA. Also runs followup script `format_cellbender_output.R` to output sparse matrices from h5.
3. `run_scrublet.sh` runs python script `run_scrublet.py` to detect doublets.


## snATAC-seq
``scripts/step1_atac_archr.sh`` does the following:
1. Downloads the fragment files from this [cart](https://www.encodeproject.org/carts/enc4_mouse_snatac/) for all 5 tissues using bash syntax.
2. Calls ``scripts/make_archr_proj.R`` which makes ArchR Arrow files, initializes ArchR Project, adds metadata, and filters cells. Cells must be present in filtered RNA object as well as passing ArchR filters (minTSS = 4, minFrags = 1000, filterDoublets).
4. Continue ArchR processing by dimensionality reduction, clustering, and plotting resulting UMAPs.

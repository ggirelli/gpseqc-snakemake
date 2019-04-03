Description of folder structure and file naming in `centrality_by_seq`, which contains the centrality scores estimated by GPSeq experiments (i.e., sequencing data analysis).

### Folder structure

* Seq: whole genome binning
    + BXXX[.exp]: sequencing run ID
        - condition(s): masking performed on the output (standard)
            + ...
        - condition(s)_noMask: no masking on the output
            + ...
* FISH: centrality estimated around probes
    + BXXX[.exp]: sequencing run ID
        - condition(s): masking performed on the output (standard)
            + ...
        - condition(s)_noMask: no masking on the output
            + ...

### Folder names

Conditions are specified in the folder names in a similar way as when specifying pages to be printed. For example:

* 1-15 means all conditions between 1 and 15 minutes (i.e., 1', 5', 10' and 15').
* 1,15 means 1' and 15' conditions only

### Filenames

Information on binning size, step, grouping and outlier removal is in the single filenames. The size of binning around probes is also specified in the single filenames (e.g., 1Mbprobe, 10Mbprobe,...).

* `B72.control.estimated.bins.size1000000.step100000.group10000.csm3.rmOutliers_chi2.rmAllOutliers.tsv`
* `B48.1Mbprobe.ranked.customBins.group10000.csm3.rmOutliers_chi2.rmAllOutliers.tsv`

If the experiments have a control (i.e., LamB_KD,...) it is coded in the filename. For example if "B10" has both a "Control" and a "LamB_KD", the files should start with "B10.ctrl." and "B10.LamB_KD.", respectively.

## Supplementary files/folders

* **snakemake**: snakemake configfile and snakefile for centrality estimation.
* **mkRecapTable.R**: Rscript to generate a recap table (to be stored on the common GoogleDrive). It produces useful warnings and errors if anything looks weird with the folder structure (as explained above).
* **cell_line_metadata.tsv**: metadata table to inform the mkRecapTable.R script (with -m option). Provide info on cell type and synchronization status.
* **GPSeq_centrality_track.tsv**: table generated by mkRecapTable.R. To be stored on the common GoogleDrive: https://docs.google.com/spreadsheets/d/12AynLiG-lleKFCJu5ov3kX97nGY6i9qmAUchSs2tYtM/edit?usp=drive_web&ouid=101270527627132364088
* **README.txt**: this file.
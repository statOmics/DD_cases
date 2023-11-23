# Differential detection project: case study

This repository contains the code associated with the case studies from our paper "Differential detection workflows for multi-sample single-cell RNA-seq data" (TODO: add link). 

## COVID case study

The COVID dataset from [Stephenson *et al.* (2021)](https://doi.org/10.1038/s41591-021-01329-2) can be downloaded programatically by running the first script of the Covid case study, 1_CovidCase_getdata.Rmd. 

Alternatively, if you would like to avoid downloading these data, you can download intermediate data objects from Zenodo (TODO: add link) and place these in the `Covid/objects/` directory. The file `sce_Covid_Bcells.rds` allows you to start the analysis from script 2_CovidCase_DGE.Rmd. Downloading the four other files from Zenodo and placing these in the `objects` directory allows you to start the analysis from scrip 4_CovidCase_downstream.Rmd.

## Lupus case study

The Lupus data from [Perez *et al.* (2021)](https://doi.org/10.1126/science.abf1970) can be
downloaded from GEO using accession number
[GSE174188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188). To run the case study code from start to end, download the `GSE174188_CLUES1_adjusted.h5ad.gz ` file and place it in the `Lupus/data-raw/` directory. Next, run the scripts 00-anndata-to-SCE.R, 01_lupus_prepare.Rmd and 02_lupus_case.Rmd in order. 

Alternatively, if you would like to avoid downloading `GSE174188_CLUES1_adjusted.h5ad.gz ` (11.7Gb), you can download intermediate data objects from Zenodo (TODO: add link) and place these in the `Lupus/objects/` directory. The file 'lupus-SCE-cleaned.rds' allows you to start the analysis from script 01_lupus_prepare.Rmd, the file `sce_sep_filt.rds` allows you to start the analysis from script 02_lupus_case.Rmd. 


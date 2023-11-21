# Differential detection project: case study

This repository contains the code associated with the case studies from our paper "Differential detection workflows for multi-sample single-cell RNA-seq data" (TODO: add link). 

## Downloading raw data

### Lupus data

The Lupus data from [Perez *et al.* (2021)](https://doi.org/10.1126/science.abf1970) can be
downloaded from GEO using accession number
[GSE174188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188).

Download the
`GSE174188_CLUES1_adjusted.h5ad.gz ` file and place it in the `data-raw/` directory under `benchmarks/lupus/`.

### COVID data

The COVID dataset from [Stephenson *et al.* (2021)](https://doi.org/10.1038/s41591-021-01329-2) should not be 
downloaded manually, as it will be downloaded automatically when running the COVID benchmarks. 

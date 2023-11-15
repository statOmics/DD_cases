## Convert lupus data to SCE and clean metadata
## ============================================

#renv::restore("../..")

## Directory setup
here_root <- "benchmarks/lupus"
here::i_am(file.path(here_root, "scripts", "00-anndata-to-SCE.R"))

in_dir <- here::here(here_root, "data-raw")
stopifnot(dir.exists(in_dir))
fs::dir_create(out_dir <- here::here(here_root, "data"))

out_file <- file.path(out_dir, "lupus-SCE-cleaned.rds")

## Check zellkonverter version
suppressPackageStartupMessages({
    library(zellkonverter)
    library(scuttle)
    library(dplyr)
})

(zell_version <- packageVersion("zellkonverter"))
stopifnot(numeric_version(zell_version) >= numeric_version("1.11.1"))

print(packageVersion("dplyr"))
print(packageVersion("scuttle"))
print(packageVersion("zellkonverter"))

## Load raw data as HDF5-backed
lupus_file <- file.path(in_dir, "GSE174188_CLUES1_adjusted.h5ad")
stopifnot(file.exists(lupus_file))

lupus <- readH5AD(lupus_file, use_hdf5 = TRUE, raw = TRUE)
stopifnot(is(assay(lupus, "X"), "DelayedMatrix"))


# Make SingleCellExperiment object ----------------------------------------

raw <- altExp(lupus, "raw")
sce <- SingleCellExperiment(
    assays = list(counts = assay(raw)),
    rowData = rowData(raw),
    colData = colData(lupus),
    metadata = metadata(lupus),
    reducedDims = list(
        PCA = reducedDim(lupus, "X_pca"),
        UMAP = reducedDim(lupus, "X_umap")
    )
)
names(rowData(sce))[2] <- "feature_types"

## Add 'sample' column: combination of 'batch_cov' and 'ind_cov'
sce$sample <- forcats::fct_cross(sce$batch_cov, sce$ind_cov, sep = ":")

## Use more descriptive name for pop_cov
sce$Ethnicity <- sce$pop_cov

## Replace missing `ct_cov` entries with `cg_cov`
sce$orig.ct_cov <- sce$ct_cov
no_ct <- is.na(sce$ct_cov)
sce$ct_cov <- as.character(sce$ct_cov)
colData(sce)[no_ct, "ct_cov"] <- as.character(colData(sce)[no_ct, "cg_cov"])
sce$ct_cov <- as.factor(sce$ct_cov)


# Sample-level metadata ---------------------------------------------------

## Create sample-level metadata table
sample_metadata <- as.data.frame(colData(sce)) %>%
    dplyr::distinct(sample, .keep_all = TRUE) %>%
    ## Keep only relevant columns
    select(
        sample, batch_cov, ind_cov, Processing_Cohort,
        Status, SLE_status,
        Age, Sex, Ethnicity,
    )
rownames(sample_metadata) <- sample_metadata$sample

## Add to metadata(lupus)
metadata(sce) <- list(sample_metadata = sample_metadata)

## add rownames
rownames(sce) <- rowData(sce)$gene_ids

# Export data -------------------------------------------------------------

## Save data as HDF5-backed SCE This will save the object with pointers to the
## HDF5 file from which the object is derived.
## See also https://github.com/bioconductor/tenxbraindata#saving-computations

saveRDS(sce, out_file)

## "Reading this file back in with readRDS should work as long as the HDF5 file
## remains in its original location"

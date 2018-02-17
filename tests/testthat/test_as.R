
stopifnot(require(Seurat))

pbmc_raw <- read.table(
    file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
    as.is = TRUE
)
pbmc_small <- CreateSeuratObject(raw.data = pbmc_raw)

test_that("minimal seurat objects can be converted", {

    sce <- as(pbmc_small, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})


pbmc_small <- NormalizeData(object = pbmc_small)

test_that("normalised seurat objects can be converted", {

    sce <- as(pbmc_small, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})

pbmc_small <- ScaleData(object = pbmc_small, vars.to.regress = c("nUMI"))

test_that("scaled seurat objects can be converted", {

    sce <- as(pbmc_small, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})

pbmc_small <- FindVariableGenes(object = pbmc_small, do.plot = FALSE)

test_that("seurat objects with variable gene information can be converted", {

    sce <- as(pbmc_small, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})

pbmc_small <- RunPCA(pbmc_small)
pbmc_small <-
    RunTSNE(pbmc_small, reduction.use = "pca", dims.use = 1:5, perplexity=10)

test_that("seurat objects with reduced dimensions can be converted", {

    sce <- as(pbmc_small, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})

pbmc1 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[1:40])
pbmc2 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[41:80])
rm(pbmc_small) # clean workspace
pbmc1@meta.data$group <- "group1"
pbmc2@meta.data$group <- "group2"
pbmc_cca <- RunCCA(pbmc1,pbmc2)

test_that("seurat objects with CCA dimensions can be converted", {

    sce <- as(pbmc_cca, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})

pbmc_cca <- AlignSubspace(
    pbmc_cca,reduction.type = "cca", grouping.var = "group", dims.align = 1:2)

test_that("seurat objects with aligned CCA dimensions can be converted", {

    sce <- as(pbmc_cca, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})

pbmc_cca <- RunPCA(pbmc_cca)

test_that("seurat objects with aligned CCA and PCA can be converted", {

    sce <- as(pbmc_cca, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})

pbmc_cca <- RunTSNE(
    pbmc_cca, reduction.use = "pca", dims.use = 1:5, perplexity=10)

test_that("seurat objects with aligned CCA and PCA can be converted", {

    sce <- as(pbmc_cca, "SingleCellExperiment")

    expect_s4_class(sce, "SingleCellExperiment")

})

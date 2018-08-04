#' Convert \code{seurat} to \code{SingleCellExperiment} objects
#'
#' @name as
#'
#' @param from A \code{\linkS4class{seurat}} object.
#'
#' @importFrom methods slot
#'
#' @return A \code{\linkS4class{SingleCellExperiment}} object.
#'
#' @examples
#' # Example data ----
#'
#' stopifnot(require(Seurat))
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
#'   as.is = TRUE
#' )
#' pbmc_small <- CreateSeuratObject(raw.data = pbmc_raw)
#' pbmc_small
#'
#' # Example ----
#'
#' as(pbmc_small, "SingleCellExperiment")
#'
setAs("seurat", "SingleCellExperiment", function(from)
{
    makeSingleCellExperimentFromSeurat(from)
})

# don't export
makeSingleCellExperimentFromSeurat <- function(from){
    # @raw.data is an absolute requirement: it is expected
    # Warning: there can be a dimension mismatch between raw.data and
    # subsequent filtered assay-related slots. This will need filtering
    seurat_raw_data <- slot(from, "raw.data")

    # @data is initialisedby CreateSeuratObject(): it is expected
    seurat_data <- slot(from, "data")

    # @cells.use contains the identifier of the cells that passed QC
    # in CreateSeuratObject()
    seurat_cells_use <- slot(from, "cell.names")

    # Seurat does not store the names of genes that passed QC
    # it directly subsets the @raw.data and @data slots
    seurat_genes_use <- rownames(slot(from, "data"))

    # The @ident slot is a useful one to keep as @colData
    seurat_ident <- slot(from, "ident")

    # ident, nGene, and nUMI in @meta.data are good @colData
    # take other metadata too, but make sure they appear after
    expected_metadata <- c("nGene","nUMI","orig.ident")
    seurat_metadata <- slot(from, "meta.data")
    metadata_seurat_order <- match(expected_metadata, colnames(seurat_metadata))
    metadata_other_order <- which(!colnames(seurat_metadata) %in% c("nGene","nUMI","orig.ident"))
    metadata_all_order <- c(metadata_seurat_order, metadata_other_order)

    # Make a SingleCellExperiment object from the expected minimal data above
    sce <- SingleCellExperiment(
        assays = list(
            raw.data = seurat_raw_data[seurat_genes_use, seurat_cells_use],
            data = seurat_data[seurat_genes_use, seurat_cells_use]
        ),
        colData = DataFrame(
            seurat_metadata[,metadata_all_order]
        )
    )

    # @scale.data may be NULL, in which case the following line does nothing
    assay(sce, "scale.data") <- slot(from, "scale.data")[seurat_genes_use, seurat_cells_use]

    # If computed, initialise @rowData with the variable genes
    if (length(slot(from, "var.genes")) > 0){
        rowData(sce) <- DataFrame(
            var.genes = rownames(sce) %in% slot(from, "var.genes"),
            hvg.info = slot(from, "hvg.info")[rownames(sce),]
        )
    }

    # Add available reduce dimension coordinates
    for (dr_name in names(slot(from, "dr"))){
        reducedDim(sce, dr_name) <- slot(slot(from, "dr")[[dr_name]], "cell.embeddings")
    }

    return(sce)
}



### Script to cluster and annotate scATAC-seq data with Signac 
# need to install MACS2 as well: pip install MACS2 --user



### User defined variables

# unpack variables passed from parent shell script
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

for (e in args) {
  argname <- e[1]
  argval <- e[2]
  # regular expression to delete initial \" and trailing \"
  argval <- gsub("(^\\\"|\\\"$)", "", argval)
  assign(argname, argval)
}



### packages

if (!require("Signac")) install.packages("Signac")
if (!require("Seurat")) install.packages("Seurat")
if (!require("Matrix")) install.packages("Matrix")
if (!require("patchwork")) install.packages("patchwork")
if (!require("GenomicRanges")) install.packages("GenomicRanges")
if (!require("rtracklayer")) BiocManager::install("rtracklayer")
if (!require("EnsDb.Hsapiens.v86")) BiocManager::install("EnsDb.Hsapiens.v86")
if (!require("EnsDb.Mmusculus.v79")) BiocManager::install("EnsDb.Mmusculus.v79")
if (!require(harmony)) devtools::install_github("immunogenomics/harmony")

if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("viridis")) install.packages("viridis")
if (!require("ggplot2")) BiocManager::install("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")

if (!require("remotes")) install.packages("remotes")
if (!require("SeuratDisk")) remotes::install_github("mojaveazure/seurat-disk")
if (!require("SeuratObject")) install.packages("SeuratObject")
if (!require("sp")) install.packages("sp")

if (!require("MuDataSeurat")) remotes::install_github("pmbio/MuDataSeurat")

set.seed(0)




### functions

## function to convert between AnnData and Seurat objects
## converting from h5ad only works if anndata<0.8.0 has been used to create h5ad!!!
# obj_file: AnnData or Seurat object filepath (AnnData .h5ad file or Seurat RDS file)
# data_name: Name of Seurat object in rdata workspace, 
# 	      only used if obj_file ends with ".rdata"
# destination: Indicates if object should be converted to AnnData or Seurat format
# out_dir: Directory that the converted object will be written to
# out_name: Name suffix for the converted object
convert_between_seurat_anndata <- function(obj_file, data_name=NULL, 
                                           destination=c("anndata","seurat"), 
                                           out_dir, out_name) {
  
  if (destination=="anndata") {
    
    if (file.exists(paste(out_dir, "/", out_name, ".h5ad", sep=""))) {
      print(paste(out_dir, "/", out_name, ".h5ad already exists, not converting", sep=""))
    } else {
      print(paste("Converting", obj_file, "to AnnData"))
      
      # read in the object
      if (unlist(strsplit(obj_file, ".", fixed=TRUE))[-1] == "rdata") {
        load(obj_file)
        obj <- eval(as.name(data_name))
      } else {
        obj <- readRDS(obj_file)
      }
      print(obj)
      
      SaveH5Seurat(obj, filename = paste(out_dir, "/", out_name, ".h5Seurat", sep=""))
      Convert(paste(out_dir, "/", out_name, ".h5Seurat", sep=""), dest = "h5ad")
      file.remove(paste(out_dir, "/", out_name, ".h5Seurat", sep=""))
      
      return(obj)
    }
  
  } else if (destination=="seurat") {
    if (file.exists(paste(out_dir, "/", out_name, sep=""))) {
      print(paste(out_dir, "/", out_name, " already exists, not converting", sep=""))
    } else {
      print(paste("Converting", obj_file, "to Seurat"))
      
      # extract input file directory and name
      filesplit <- unlist(strsplit(obj_file, "/"))
      in_dir <- paste(filesplit[1:(length(filesplit)-1)], collapse="/")
      filename <- filesplit[length(filesplit)]
      filenamesplit <- unlist(strsplit(filename, "\\."))
      filename_wo_suffix <- paste(filenamesplit[1:(length(filenamesplit)-1)], collapse="_")
      
      Convert(obj_file, dest = "h5seurat")
      obj <- LoadH5Seurat(paste(in_dir, "/", filename_wo_suffix, ".h5seurat", sep=""))
      
      # save as RDS file
      saveRDS(obj, paste(out_dir, "/", out_name, sep=""))
      file.remove(paste(in_dir, "/", filename_wo_suffix, ".h5seurat", sep=""))
      
      return(obj)
    }
  }
}




## subfunction to create peak count matrix from fragment file
create_peak_count_matrix <- function(fragment_file_path, frag_cutoff=1000, features=NULL,
                                     species=c("human","mouse","zebrafish"), out_dir, out_name) {
  
  # Define cells
  # If you already have a list of cell barcodes to use you can skip this step
  total_counts <- CountFragments(fragment_file_path)
  # Change frag_cutoff depending on your dataset!
  barcodes <- total_counts[total_counts$frequency_count > frag_cutoff, ]$CB
  
  # Create a fragment object
  frags <- CreateFragmentObject(path = fragment_file_path, cells = barcodes)
  
  if (is.null(features)) {
    # First call peaks on the dataset
    # If you already have a set of peaks you can skip this step
    if (species=="human") {
      peaks <- CallPeaks(frags, verbose=TRUE, outdir=out_dir, effective.genome.size = 3e+09)
    } else if (species=="mouse") {
      peaks <- CallPeaks(frags, verbose=TRUE, outdir=out_dir, effective.genome.size = 2.7e+09)
    } else if (species=="zebrafish") {
      peaks <- CallPeaks(frags, verbose=TRUE, outdir=out_dir, effective.genome.size = 1.5e+09)
    } else {
      print("Could not call peaks, set species as 'human', 'mouse' or 'zebrafish'")
    }
    
    # subset to main chromosomes
    peaks <- peaks[seqnames(peaks) %in% standardChromosomes(peaks)]
    
    # Quantify fragments in each peak
    counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
  
  } else {
    # Quantify fragments in each peak
    counts <- FeatureMatrix(fragments = frags, features = features, cells = barcodes)
  }
  
  saveRDS(counts, paste(out_dir, "/", out_name, "_signac_feature_matrix.rds", sep=""))

  return(counts)
}



## function to create signac object from single fragment file
create_signac <- function(fragment_file_path, counts, species=c("human","mouse","zebrafish"), 
                          metadata_file_path=NULL, out_name, out_dir) {
  
  # adjust chromosome names in feature matrix
  if (grepl("chr", rownames(counts)[1], fixed = TRUE)) {
    orig_chrom_style <- "ucsc"
  } else {
    orig_chrom_style <- "ensembl"
  }

  if (orig_chrom_style=="ensembl") {
    rownames(counts) <- paste("chr", rownames(counts), sep="")
    rownames(counts) <- gsub("chrMT", "chrM", rownames(counts))
  }

  # create signac object
  chrom_assay <- CreateChromatinAssay(counts = counts,
                                      sep = c(":", "-"),
                                      fragments = fragment_file_path,
                                      min.cells = 10,
                                      min.features = 200)

  if (!is.null(metadata_file_path)) {
    metadata <- read.table(file = metadata_file_path,
                           header = TRUE,
                           sep = "\t",
                           row.names = 1)
  
    signac_obj <- CreateSeuratObject(counts = chrom_assay,
                                     assay = "ATAC",
                                     meta.data = metadata)

    # subset object to cells contained in metadata file
    signac_obj <- signac_obj[, rownames(signac_obj[[]]) %in% rownames(metadata)]

  } else {
    signac_obj <- CreateSeuratObject(counts = chrom_assay,
                                     assay = "ATAC")
  }

  # calculate total number of fragments in file
  total_fragments <- CountFragments(fragment_file_path)
  rownames(total_fragments) <- total_fragments$CB
  signac_obj$fragments <- total_fragments[colnames(signac_obj), "frequency_count"]  

  signac_obj$dataset <- out_name

  # extract gene annotations from EnsDb
  if (species=="human") {
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
 
    # change to UCSC style
    seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
    genome(annotations) <- "hg38"

  } else if (species=="mouse") {
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  
    # change to UCSC style
    seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
    genome(annotations) <- "mm10"

  } else if (species=="zebrafish") {
    txdb <- makeTxDbFromBiomart(dataset="drerio_gene_ensembl")
  
    # change to UCSC style
    newSeqNames <- paste0('chr', seqlevels(txdb))
    names(newSeqNames) <- seqlevels(txdb)
    txdb <- renameSeqlevels(txdb, newSeqNames)

    # generate GRanges list object, listing exons per gene
    tmp <- exonsBy(txdb, by="gene")

    # combine overlapping exons for individual genes
    annotations <- reduce(tmp)
  }
  
  # add the gene information to the object
  Annotation(signac_obj) <- annotations

  print(signac_obj)
  print(granges(signac_obj))
  
  # export peak locations
  peaks_df = as(granges(signac_obj), "data.frame")
  write.table(peaks_df, paste(out_dir, "/", out_name, "_peaks.bed", sep=""),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  # extract cell barcodes
  cells <- Cells(signac_obj)
  write.csv(cells, paste(out_dir, "/", out_name, "_cell_barcodes.csv", sep=""),
            quote=FALSE, row.names=FALSE)
  
  # save object
  saveRDS(signac_obj, paste(out_dir, "/", out_name, "_signac.rds", sep=""))
  
  return(signac_obj)
}



## function to create Signac object from multiple input fragment files
# input_file_list: list of fragment file filepaths
# sample_list: list of fragment file identifiers
create_consensus_signac <- function(input_file_list, sample_list,
                                    species=c("human","mouse","zebrafish"),  
                                    metadata_file_path=NULL,
                                    out_dir, out_name) {
  
  # create raw Signac objects with individual peak sets
  for (i in seq(1,length(input_file_list))) {

    if (!file.exists(paste(out_dir, "/", sample_list[i], "_signac", sep=""))) {
      counts <- create_peak_count_matrix(fragment_file_path=input_file_list[i], frag_cutoff=1000,
                                         species=species, out_dir=out_dir, out_name=sample_list[i])

      signac_obj <- create_signac(fragment_file_path=input_file_list[i], counts=counts, species=species, 
                                  metadata_file_path=NULL,
                                  out_dir=out_dir, out_name=sample_list[i])
    }
  
    # read in peak sets + convert to genomic ranges
    peaks <- makeGRangesFromDataFrame(read.table(file = paste(out_dir, "/", sample_list[i], 
                                                              "_peaks.bed", sep=""), header=TRUE))

    # subset to main chromosomes
    peaks <- peaks[seqnames(peaks) %in% standardChromosomes(peaks)]
    
    # combine peaks into single granges object
    if (i==1) {
      all_peaks_combined <- peaks
    } else {
      all_peaks_combined <- c(all_peaks_combined, peaks)
    }
  }
  print(paste("Peaks before reducing:", nrow(as.data.frame(all_peaks_combined))))

  # fuse overlapping peaks
  combined.peaks <- reduce(x = all_peaks_combined)
  print(paste("Peaks after reducing:", nrow(as.data.frame(combined.peaks))))
  
  # Filter out bad peaks based on length
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
  
  for (i in seq(1,length(input_file_list))) {
    if (!file.exists(paste(out_dir, "/", sample_list[i], "_cons_signac", sep=""))) {
      print(paste("Building consensus peak Signac object for ", sample_list[i], sep=""))
      counts <- create_peak_count_matrix(fragment_file_path=input_file_list[i], frag_cutoff=1000,
                                         features=combined.peaks, species=species,
                                         out_dir=out_dir, out_name=paste(sample_list[i],"_cons",sep=""))
 
      # adjust chromosome names in feature matrix
      if (grepl("chr", rownames(counts)[1], fixed = TRUE)) {
        orig_chrom_style <- "ucsc"
      } else {
        orig_chrom_style <- "ensembl"
      }

      if (orig_chrom_style=="ensembl") {
        rownames(counts) <- paste("chr", rownames(counts), sep="")
        rownames(counts) <- gsub("chrMT", "chrM", rownames(counts))
      }
      
      signac_obj <- create_signac(fragment_file_path=input_file_list[i], counts=counts, species=species, 
                                  metadata_file_path=metadata_file_path,
                                  out_dir=out_dir, out_name=paste(sample_list[i],"_cons",sep=""))

      # subset object to cells contained in metadata file
      signac_obj <- signac_obj[, rownames(signac_obj[[]]) %in% rownames(metadata)]

    } else {
      print(paste("Using existing consensus peak Signac object for ", sample_list[i], sep=""))
      signac_obj <- readRDS(paste(out_dir, "/", sample_list[i], "_cons_signac", sep=""))
    }
    
    signac_obj$sample_id <- sample_list[i]
    signac_obj$dataset <- out_name
    
    # calculate total number of fragments in file
    total_fragments <- CountFragments(input_file_list[i])
    rownames(total_fragments) <- total_fragments$CB
    signac_obj$fragments <- total_fragments[colnames(signac_obj), "frequency_count"]

    print(head(signac_obj[[]]))
    print(granges(signac_obj[[]]))    
    print(length(Cells(signac_obj)))

    # merge new signac object into combined object
    if (i==1) {
      combined_signac_obj <- signac_obj
    } else {
      tmp <- merge(x = combined_signac_obj, y = signac_obj)
      combined_signac_obj <- tmp
    }
  }

  # save object
  saveRDS(combined_signac_obj, paste(out_dir, "/", out_name, "_signac_combined.rds", sep=""))
  
  print(combined_signac_obj)
  
  return(combined_signac_obj)
}



## function to do QC on signac object
qc_signac <- function(signac_obj, species=c("human","mouse","zebrafish"), 
                      min_count=3000, max_count=30000, out_dir, out_name) {
  
  # compute nucleosome signal score per cell
  signac_obj <- NucleosomeSignal(object = signac_obj)
  
  # compute TSS enrichment score per cell
  signac_obj <- TSSEnrichment(object = signac_obj, assay = "ATAC", fast = FALSE)
  
  # add fraction of reads in peaks
  signac_obj <- FRiP(object = signac_obj, assay = "ATAC", total.fragments = "fragments")
  signac_obj$pct_reads_in_peaks <- signac_obj$FRiP * 100
  
  # add blacklist ratio
  if (species == "human") {
    signac_obj$blacklist_ratio <- FractionCountsInRegion(object = signac_obj, assay = 'ATAC',
                                                         regions = blacklist_hg38)
  } else if (species == "mouse") {
    signac_obj$blacklist_ratio <- FractionCountsInRegion(object = signac_obj, assay = 'ATAC',
                                                         regions = blacklist_mm10)
  }
  
  # plot density of peaks over TSS enrichment
  #gg <- DensityScatter(signac_obj, x = 'nCount_peaks', y = 'TSS.enrichment', 
  #                     log_x = TRUE, quantiles = TRUE)
  #ggsave(gg, filename=paste(out_dir, "/", out_name, "_density_scatter.pdf", sep=""), device="pdf", 
  #       dpi=300, useDingbats=FALSE)
  
  # plot TSS enrichment
  signac_obj$high.tss <- ifelse(signac_obj$TSS.enrichment > 3, 'High', 'Low')
  gg <- TSSPlot(signac_obj, group.by = 'high.tss') + NoLegend()
  ggsave(gg, filename=paste(out_dir, "/", out_name, "_tss_enrichment.pdf", sep=""), device="pdf", 
         dpi=300, useDingbats=FALSE)
  
  # plot fragment size distribution
  signac_obj$nucleosome_group <- ifelse(signac_obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  #gg <- FragmentHistogram(object = signac_obj, group.by = 'nucleosome_group')
  #ggsave(gg, filename=paste(out_dir, "/", out_name, "_fragment_size.pdf", sep=""), device="pdf", 
  #       dpi=300, useDingbats=FALSE)
  
  # plot QC stats
  gg <- VlnPlot(object = signac_obj, features = c('nCount_ATAC', 'TSS.enrichment', 
                                                  'blacklist_ratio', 'nucleosome_signal', 
                                                  'pct_reads_in_peaks'), pt.size = 0, ncol = 5)
  ggsave(gg, filename=paste(out_dir, "/", out_name, "_qc_violin.pdf", sep=""), device="pdf", 
         dpi=300, useDingbats=FALSE, width=14)
  
  # exclude bad cells and peaks
  if (species %in% c("human","mouse")) {
    signac_obj <- subset(x = signac_obj, subset = nCount_ATAC > min_count &
                                                  nCount_ATAC < max_count &
                                                  pct_reads_in_peaks > 15 &
                                                  blacklist_ratio < 0.05 &
                                                  nucleosome_signal < 4 &
                                                  TSS.enrichment > 3)
  } else {
    signac_obj <- subset(x = signac_obj, subset = nCount_ATAC > min_count &
                                                  nCount_ATAC < max_count &
                                                  pct_reads_in_peaks > 15 &
                                                  nucleosome_signal < 4 &
                                                  TSS.enrichment > 3)
  }
  
  # save object
  saveRDS(signac_obj, paste(out_dir, "/", out_name, "_signac_qc.rds", sep=""))
  
  print(signac_obj)
  
  return(signac_obj)
}



## function to cluster signac object
cluster_signac <- function(signac_obj, seurat_obj=NULL, integration_key=NULL, resolution=1,
                           umap_colors_rna=NULL, umap_colors_atac=NULL,
                           out_dir, out_name) {

  set.seed(0)

  signac_obj <- RunTFIDF(signac_obj)
  signac_obj <- FindTopFeatures(signac_obj, min.cutoff = 'q5')
  signac_obj <- RunSVD(signac_obj)
  
  gg <- DepthCor(signac_obj)
  ggsave(gg, filename=paste(out_dir, "/", out_name, "_component_depth.pdf", sep=""), device="pdf", 
         dpi=300, useDingbats=FALSE)
   
  # integrate dataset on integration_key
  if (!is.null(integration_key)) {
    
    print(paste("Performing harmonisation on variable:", integration_key))
    print(signac_obj)
    print(head(signac_obj[[]]))
   
    pdf(file = paste(out_dir, "/", out_name, "_harmony_run.pdf", sep=""))
      options(repr.plot.height = 2.5, repr.plot.width = 6)
      signac_obj <- RunHarmony(signac_obj, integration_key, 
                               project.dim = FALSE, 
                               assay.use = "ATAC", plot_convergence = TRUE, max_iter = 25,
                               reduction.use = "lsi")
    dev.off()
    
    # compare PCs before and after harmonisation
    pdf(file = paste(out_dir, "/", out_name, "_lsi_without_harmony.pdf", sep=""))
      options(repr.plot.height = 5, repr.plot.width = 12)
      p1 <- DimPlot(object = signac_obj, reduction = "lsi", pt.size = .1, group.by = integration_key)
      p2 <- VlnPlot(object = signac_obj, features = "LSI_2", group.by = integration_key, pt.size = .1)
      print(cowplot::plot_grid(p1,p2))
    dev.off()
    
    pdf(file = paste(out_dir, "/", out_name, "_lsi_with_harmony.pdf", sep=""))
      options(repr.plot.height = 5, repr.plot.width = 12)
      p1 <- DimPlot(object = signac_obj, reduction = "harmony", pt.size = .1, group.by = integration_key)
      p2 <- VlnPlot(object = signac_obj, features = "harmony_1", group.by = integration_key, pt.size = .1)
      print(cowplot::plot_grid(p1,p2))
    dev.off()
    
    signac_obj <- FindNeighbors(signac_obj, reduction = "harmony", dims = 2:30)
    signac_obj <- FindClusters(signac_obj, resolution = resolution, verbose = FALSE, algorithm = 3)

    # compare UMAP plots with and without integration
    signac_obj <- RunUMAP(signac_obj, reduction = 'lsi', dims = 2:40)

    gg <- DimPlot(object = signac_obj, group.by = integration_key, label = FALSE, pt.size=0.7, raster=FALSE)
    ggsave(gg, filename=paste(out_dir, "/", out_name, "_umap_atac_", integration_key, "_unintegrated.png", sep=""), 
           device="png", dpi=300, width=10)

    signac_obj <- RunUMAP(signac_obj, reduction = "harmony", dims = 2:30)
    
    gg <- DimPlot(object = signac_obj, group.by = integration_key, label = FALSE, pt.size=0.7, raster=FALSE)
    ggsave(gg, filename=paste(out_dir, "/", out_name, "_umap_atac_", integration_key, "_integrated.png", sep=""), 
           device="png", dpi=300, width=10)
    
  } else {
    signac_obj <- FindNeighbors(signac_obj, reduction = 'lsi', dims = 2:40)
    signac_obj <- FindClusters(signac_obj, resolution = resolution, verbose = FALSE, algorithm = 3)

    signac_obj <- RunUMAP(signac_obj, reduction = 'lsi', dims = 2:40)

    gg <- DimPlot(object = signac_obj, group.by = integration_key, label = FALSE, pt.size=0.7, raster=FALSE)
    ggsave(gg, filename=paste(out_dir, "/", out_name, "_umap_atac_", integration_key, "_unintegrated.png", sep=""), 
           device="png", dpi=300, width=10)
  }
  
  gg <- DimPlot(object = signac_obj, label = TRUE, pt.size=0.7, raster=FALSE) + NoLegend()
  ggsave(gg, filename=paste(out_dir, "/", out_name, "_umap_atac_clusters.png", sep=""), device="png", 
         dpi=300)
  
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  gene.activities <- GeneActivity(signac_obj)
  
  signac_obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
  signac_obj <- NormalizeData(object = signac_obj, assay = 'RNA',
                              normalization.method = 'LogNormalize',
                              scale.factor = median(signac_obj$nCount_RNA))
  print(head(rownames(signac_obj)))
  print(head(rownames(seurat_obj)))
  
  # transfer cell labels
  if (!is.null(seurat_obj)) {
    transfer.anchors <- FindTransferAnchors(reference = seurat_obj, query = signac_obj, 
                                            query.assay = "RNA", reduction = 'cca',
                                            features = rownames(seurat_obj))
    predicted.labels <- TransferData(anchorset = transfer.anchors, 
                                     refdata = seurat_obj$cell_type,
                                     weight.reduction = signac_obj[['lsi']],
                                     dims = 2:40)
  
    signac_obj <- AddMetaData(object = signac_obj, metadata = predicted.labels)
  
    # plot umap of cell type labels
    plot1 <- DimPlot(object = seurat_obj, group.by = 'cell_type', pt.size=0.7, cols=umap_colors_rna,
                     label = FALSE, repel = TRUE, raster=FALSE) + ggtitle(paste("scRNA-seq", out_name))
    ggsave(plot1, filename=paste(out_dir, "/", out_name, "_umap_celltype_rna.png", sep=""), device="png", 
           dpi=300, width=10)
  
    plot2 <- DimPlot(object = signac_obj, group.by = 'predicted.id', pt.size=0.7, cols=umap_colors_atac,
                     label = FALSE, repel = TRUE, raster=FALSE) + ggtitle(paste("scATAC-seq", out_name))
    ggsave(plot2, filename=paste(out_dir, "/", out_name, "_umap_celltype_atac.png", sep=""), device="png", 
           dpi=300, width=10)
  
    # replace each label with its most likely prediction
    for(i in levels(signac_obj)) {
      cells_to_reid <- WhichCells(signac_obj, idents = i)
      newid <- names(which.max(table(signac_obj$predicted.id[cells_to_reid])))
      Idents(signac_obj, cells = cells_to_reid) <- newid
    }
    signac_obj$cell_type <- Idents(signac_obj)
  
    gg <- DimPlot(object = signac_obj, pt.size=0.7, cols=umap_colors_atac,
                  label = FALSE, repel = TRUE, raster=FALSE) + ggtitle(paste("scATAC-seq", out_name))
    ggsave(gg, filename=paste(out_dir, "/", out_name, "_umap_celltype_atac_reid.png", sep=""), device="png", 
           dpi=300, width=10)
  }
  
  # save object
  saveRDS(signac_obj, paste(out_dir, "/", out_name, "_signac_qc_analysed.rds", sep=""))
  
  print(signac_obj)
  
  return(signac_obj)
}



## function to plot UMAP
# idents_list: list of layers to overlay onto the UMAP clustering
# resolution: clustering resolution used, not the resolution of the saved images 
# if labels=TRUE, text labels will be displayed on UMAP and plot will be saved as pdf file
# if labels=FALSE, text labels will not be displayed and plot will be saved as png file
# file_type must be compatible with ggsave: "pdf", "png", etc.
umap_seurat <- function(seurat_obj, out_dir, project_name, idents_list, labels=TRUE, 
                        resolution=0.4, file_type="pdf") {
  
  for (ident in idents_list) {
    
    if (labels == TRUE) {
      if (ident == "seurat_clusters") {
        resolution <- resolution * 100
        filename <- paste(out_dir, "/", project_name, "_umap_res_0_", resolution, 
                          ".", file_type, sep="")
      } else {
        filename <- paste(out_dir, "/", project_name, "_umap_", ident, ".", 
                          file_type, sep="")
      }

      if (is.numeric(seurat_obj[[ident]][1,]) | is.integer(seurat_obj[[ident]][1,])) {
          gg <- FeaturePlot(seurat_obj, reduction = "umap", features=ident,
                            label = FALSE, pt.size = 0.4)
      } else {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident,
                        label = TRUE, pt.size = 0.4) + NoLegend() 
      }

      ggsave(gg, filename=filename, device=file_type, dpi=300)

    } else if (labels == FALSE) {
      if (ident == "seurat_clusters") {
        resolution <- resolution * 100
        filename <- paste(out_dir, "/", project_name, "_umap_res_0_", resolution, 
                          "_no_labels.", file_type, sep="")
      } else {
        filename <- paste(out_dir, "/", project_name, "_umap_", ident, "_no_labels.", 
                          file_type, sep="")
      }

      if (is.numeric(seurat_obj[[ident]][1,]) | is.integer(seurat_obj[[ident]][1,])) {
        gg <- FeaturePlot(seurat_obj, reduction = "umap", features=ident,
                          label = FALSE, pt.size = 0.4)
      } else {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident,
                        label = FALSE, pt.size = 0.4)     
      }

      ggsave(gg, filename=filename, device=file_type, dpi=300, width=11)
    }
  }

  return(seurat_obj)
}





## function to plot stacked bar chart of normalised cell counts
# x_axis_columns: list of one or multiple metadata categories that data will be grouped by on the x axis,
#                   category names need to be column names in metadata,
#                   columns should contain strings
# fill_column: name of metadata column that stacked bars will be fill coloured by,
#              column should contain strings
# x_axis_levels: optional list of x_axis_columns combinations to plot (separator "__"), 
#                to change order of plotting on x axis
# fill_levels: optional list of fill_column values,
#              to change order of plot fill and legend
metadata_stacked_bar <- function(metadata, x_axis_columns, fill_column,    
                                 x_axis_levels=NULL, fill_levels=NULL,
                                 out_dir, out_name) {
  
  # group and summarise metadata on columns in x_axis_columns and on fill_column
  if (length(x_axis_columns) == 1) {
    plot_df <- metadata %>% group_by(eval(as.name(x_axis_columns)), 
                                      eval(as.name(fill_column))) %>% summarise(count=n())
    normalisation_df <- metadata %>% group_by(eval(as.name(x_axis_columns))) %>% summarise(total_count=n())
  } else if (length(x_axis_columns) > 1) {
    metadata$combined_cat <- apply(metadata[,eval(x_axis_columns)], 1, paste, collapse="__")
    plot_df <- metadata %>% group_by(combined_cat, 
                                      eval(as.name(fill_column))) %>% summarise(count=n())
    normalisation_df <- metadata %>% group_by(combined_cat) %>% summarise(total_count=n())
  }

  colnames(plot_df)[1:2] <- c("cat_x", "cat_fill")
  colnames(normalisation_df)[1] <- "cat_x"

  # normalise counts per "cat_x" value
  plot_df_1 <- merge(plot_df, normalisation_df, by="cat_x")
  plot_df_1$norm_count <- plot_df_1$count / plot_df_1$total_count

  # optionally relevel x axis elements 
  if (!is.null(x_axis_levels)) {
    plot_df_1$cat_x <- factor(plot_df_1$cat_x, levels = x_axis_levels)
  }

  # optionally relevel y axis elements 
  if (!is.null(fill_levels)) {
    plot_df_1$cat_fill <- factor(plot_df_1$cat_fill, levels = fill_levels)
  }


  # create stacked bar plot
  colours <- c("#ff7f0e","#2ca02c","#d62728","#9467bd",
               "#8c564b","#e377c2","#bcbd22","#1f77b4",
               "#98df8a","#ffbb78","#aec7e8","#17becf",
               "#ff9896","#c5b0d5","#c49c94","#f7b6d2",
               "#dbdb8d","#9edae5")
  #colours <- c("dodgerblue", "darkgreen", "deeppink3", "lightgoldenrod", "saddlebrown", 
  #             "grey60", "springgreen2", "darkorange", "slateblue4", "maroon4",
  #             "black", "firebrick2", "bisque3", "darkblue", "indianred3",
  #             "turquoise2", "grey30")

  gg <- ggplot(plot_df_1, aes(x=cat_x,y=norm_count)) +
    geom_bar(aes(fill=cat_fill),stat="identity", width=0.8) +
    scale_fill_manual(values=colours) +
    theme_minimal(base_size = 22) +
    theme(
    axis.text.x = element_text(size = 22, colour = "black",face = "plain", angle=30, hjust=1),
    axis.text.y = element_text(size = 22, colour = "black",face = "plain"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "plain", vjust = 0),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "right",
    legend.text = element_text(size = 18, colour = "black"),
    legend.title = element_blank(),
    plot.margin = margin(0, 1, 0, 5, "cm"),
    strip.text.x = element_text(size=16, face="plain"),
    strip.background = element_blank()) +
    scale_y_continuous(name="Relative count", limits=c(0, 1), breaks=seq(0,1,by=.2)) +
    scale_x_discrete() +
    labs(title="")
  
  ggsave(gg, filename=paste(out_dir, "/", out_name, "_", fill_column, "_vs_",
                        paste(x_axis_columns,collapse="_"), "_stacked_bar.pdf",sep=""),
         device="pdf", dpi=300, useDingbats=FALSE, height=9, width=15)

  return(plot_df_1)
}





### Analysis

leiden_res <- as.numeric(leiden_res)



# create output directory if it does not exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}



## read in scRNA-seq object
if (file.exists(scRNA_path)) {
  if (unlist(strsplit(basename(scRNA_path), split="\\."))[-1] == "rds") {
    seurat_obj <- ReadRDS(scRNA_path)
  } else if (unlist(strsplit(basename(scRNA_path), split="\\."))[-1] == "h5ad") {
    seurat_obj <- ReadH5AD(scRNA_path)

#   out_name_rna <- paste(out_name, "_scanpy", sep="")
#   seurat_obj <- convert_between_seurat_anndata(obj_file=scRNA_path, data_name=NULL, 
#                                                destination="seurat", 
#                                                out_dir=out_dir, out_name=project_name)
  }
} else {
  seurat_obj <- NULL
}



## generate Signac object
if (!file.exists(paste(out_dir, "/", project_name, "_signac.rds", sep=""))) {
  if (!file.exists(paste(out_dir, "/", project_name, "_signac_feature_matrix.rds", sep=""))) {
    counts <- create_peak_count_matrix(fragment_file_path=in_file, frag_cutoff=1000,
                                     species=species, out_dir=out_dir, out_name=project_name)
  } else {
    counts <- readRDS(paste(out_dir, "/", project_name, "_signac_feature_matrix.rds", sep=""))
  }

  signac_obj <- create_signac(fragment_file_path=in_file, 
                              counts=counts, 
                              species=species,
                              metadata_file_path=metadata,
                              out_dir=out_dir, out_name=project_name)
} else {
  signac_obj <- readRDS(paste(out_dir, "/", project_name, "_signac.rds", sep=""))
}



## perform quality control
if (!file.exists(paste(out_dir, "/", project_name, "_signac_qc.rds", sep=""))) {
  signac_obj <- qc_signac(signac_obj=signac_obj, species=species,
                          min_count=3000, max_count=15000,
                          out_dir=out_dir, out_name=project_name)
} else {
  signac_obj <- readRDS(paste(out_dir, "/", project_name, "_signac_qc.rds", sep=""))
}



## cluster + transfer cell labels
if (!file.exists(paste(out_dir, "/", project_name, "_signac_qc_analysed.rds", sep=""))) {
  signac_obj <- cluster_signac(signac_obj=signac_obj, seurat_obj=seurat_obj, 
                               integration_key="sample_id", resolution=leiden_res,
                               umap_colors_rna=NULL, 
                               umap_colors_atac=NULL,
                               out_dir=out_dir, out_name=project_name)
} else {
  signac_obj <- readRDS(paste(out_dir, "/", project_name, "_signac_qc_analysed.rds", sep=""))
}


## plot UMAP with metadata overlays
idents_list = colnames(signac_obj[[]])
idents_list = idents_list[!idents_list %in% c("barcode", "cell_id")]
signac_obj <- umap_seurat(seurat_obj=signac_obj, out_dir=out_dir, project_name=project_name, 
                          resolution=leiden_res, idents_list=idents_list, labels=FALSE, 
                          file_type="png")


## export metadata as csv
meta_data_export <- signac_obj[[]]
write.csv(meta_data_export, paste(out_dir, "/", project_name, 
                                  "_scATACseq_annotated_metadata.csv", sep=""))


## export count matrix for pycisTopic
if (!file.exists(paste(out_dir, "/", project_name, "_counts.mtx", sep=""))) {
  writeMM(GetAssayData(object = signac_obj, slot = "counts"), 
          file=paste(out_dir, "/", project_name, "_counts.mtx", sep=""))
}



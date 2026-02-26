#!/opt/conda/bin/R
args = commandArgs(trailingOnly = TRUE)
# Parse key=value pairs
args_list = list()
i = 1
while (i <= length(args)) {
  if (startsWith(args[i], "--")) {
    key = substring(args[i], 3)           # remove leading "--"
    value = args[i + 1]                   # next argument is the value
    args_list[[key]] = value
    i = i + 2
  } else {
    i = i + 1
  }
}

print(args_list)

dataset = args_list[["batch"]]
outdir = args_list[["outdir"]]
INFILE = args_list[["INFILE"]]
doubletProportion = as.numeric(args_list[["doubletProportion"]])
numCore = as.numeric(args_list[["numCore"]])
print(INFILE)
print(doubletProportion)
print(numCore)

if (length(doubletProportion)>1) {
  stop("Error: doubletProportion must either be NULL or a single integer value!")
}

# dataset = "Eleni_Female"
# outdir = "../results/QC/"
# minSolidity = 0.5
# minArea.um2 = 20
# minnCount_RNA = 25
# doubletProportion = 0.05
# INDIR = "../data/Shared_Folder/"
# numCore = 4


findDoublets = function(data.filt, doubletProportion= NULL, numCore=1) {
    # before doing doublet detection we need to run scaling, variable gene selection and PCA, as well as UMAP for visualization. 
    # These steps will be explored in more detail in coming exercises.
    # Remove unnecessary:
    # data.filt@images = list()
    # data.filt = DietSeurat(
    #   object = data.filt,
    #   assays = "RNA",       # keep only RNA
    #   counts = TRUE,        # keep counts
    #   data = TRUE           # keep normalized data
    # )

    DefaultAssay(data.filt) = "RNA"
    data.filt = NormalizeData(data.filt)
    data.filt = ScaleData(data.filt)
    data.filt = FindVariableFeatures(data.filt)
    data.filt = RunPCA(data.filt, verbose = F, npcs = 100)

    # Then we run doubletFinder, selecting first 10 PCs and a `pK` value of 0.9. To optimize the parameters, you can run the `paramSweep` function in the package.
    # pK identification (no ground-truth)
    # pK: neighborhood size, i.e. the number of neighbors (in percentage) to consider when calculating pANN
    ## Downsample the object 
    nCellsSample = min(50000, ncol(data.filt))  # number of cells to sample
    set.seed(123)  # for reproducibility
    cells.sample = sample(Cells(data.filt), nCellsSample)
    data.sub = subset(data.filt, cells = cells.sample)
    # set seed
    set.seed(8)
    sweep.list = invisible(paramSweep(data.sub, PCs = 1:10, num.cores = numCore/2))
    sweep.stats = summarizeSweep(sweep.list)
    pdf("find.pK.plots.pdf")
      bcmvn = find.pK(sweep.stats)
    dev.off()
    # Optimal pK are selected by maximizing BCmvn (https://www.sciencedirect.com/science/article/pii/S2405471219300730)
    bcmvn.max = bcmvn[which.max(bcmvn$BCmetric),]
    optimal.pk = bcmvn.max$pK
    optimal.pk = as.numeric(levels(optimal.pk))[optimal.pk]
    print(paste0(dataset, " optimal.pk = ", optimal.pk))

    runDF = function(doubletProportion, data.filt, optimal.pk) {
        # The expected doublet rate as suggested by 10X
        # The percentage of droplets containing two or more cells is influenced by the quantity of cells loaded into the sequencing machine.
        # ref: https://uofuhealth.utah.edu/huntsman/shared-resources/gcb/htg/single-cell/genomics-10x
        cells.nr = ncol(data.filt)
        nExp.pois = round(cells.nr * doubletProportion)
        # nExp.pois

        # run DoubletFinder
        # pN: the proportion of artificial doublet cells to generate (default: 0.25)
        data.filt = doubletFinder(seu = data.filt,
                                PCs = 1:10,
                                pN = 0.10, # Artificial doublets
                                pK = optimal.pk,
                                nExp = nExp.pois)
        metadata = data.filt@meta.data
        colnames(metadata)[ncol(metadata)] = "doublet_finder"
        return(metadata)
    }

    if (doubletProportion=="null"||doubletProportion==""||is.na(doubletProportion)||is.null(doubletProportion)) {
      print("Run doubletFinder ...")
      doubletProportion_list = c(0.02, 0.05, 0.08)
      metadata = lapply(doubletProportion_list, function(d) {
        print(d)
        t = runDF(doubletProportion=d, 
              data.filt=data.filt, 
              optimal.pk=optimal.pk)
        return(t)
      }) %>% setNames(paste0("doubletProportion_", doubletProportion_list))
    } else {
      metadata = runDF(doubletProportion=doubletProportion, 
                        data.filt=data.filt, 
                        optimal.pk=optimal.pk)
    }

    return(metadata)
}
 

print("Load packages...")
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(DoubletFinder)

print(paste0("Running QC for dataset: ", dataset))
OUTDIR = paste0(outdir, "/",dataset, "/")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

dat_filtered = readRDS(INFILE)
# print("Remove potential doublets with Area.um2 & nCount_RNA that are both 8 Median Absolute Deviation (MAD) away from the median of Area.um2 & nCount_RNA")
# pdf(paste0( "./test.pdf"), h=10, w=10)
# potentialDoublets = mapply(function(f, i) {
#   median = median(dat@meta.data[[f]])
#   MAD = i * median(abs(dat@meta.data[[f]]-median))
#   potentialDoublets = dat@meta.data[[f]] > (median + MAD)
#     print(ImageDimPlot(dat, fov = dataset, cols = "red", 
#       cells = rownames(dat@meta.data)[potentialDoublets]) + 
#         ggtitle(paste0(f,">",(median + MAD))))
#   return(potentialDoublets)
# },c("nCount_RNA", "Area.um2"), c(8, 8))
# sum((potentialDoublets %>% rowSums())==2)
# print(ImageDimPlot(dat, fov = dataset, cols = "red", 
#       cells = rownames(dat@meta.data)[(potentialDoublets %>% rowSums())==2]) + 
#         ggtitle("Combined"))
# dev.off()

print("Detect doublets using doubletFinder ...")
# source("rmDoublets.R")
rm(dat)
gc()
temp = dat_filtered
temp@images = list()
temp = DietSeurat(
      object = temp,
      assays = "RNA",       # keep only RNA
      counts = TRUE,        # keep counts
      data = TRUE           # keep normalized data
)
tissue.list = SplitObject(temp, split.by = "Tissue")
metadatas = lapply(names(tissue.list), function(t) {
  print(t)
  tissue = tissue.list[[t]]
  # cells.sample = sample(Cells(tissue), 1000)
  # tissue = subset(tissue, cells = cells.sample)
  dat = findDoublets(data.filt=tissue, 
    doubletProportion= doubletProportion,
    numCore=numCore)
  return(dat)
}) %>% setNames(names(tissue.list))

print("Save doublets info ...")
paste0(OUTDIR, "/Doublets/") %>% 
    dir.create(., recursive=TRUE, showWarnings=FALSE)
saveRDS(metadatas, paste0(OUTDIR, "/Doublets/Doublet_metadatas.RDS"))

print("Reformulate the doublet inforamtion ...")
if(length(metadatas[[1]])==3) {
  metadata = lapply(names(metadatas[[1]]), function(i) {
    m = lapply(metadatas, function(tissue) {
      mdat = tissue[[i]]
      mdat[[grep("pANN_", colnames(mdat), value=TRUE)]]=NULL
      mdat
    }) %>% Reduce(rbind, .) %>% as.data.frame()
    m[[i]] = m$doublet_finder
    m[,ncol(m),drop=FALSE]
  }) %>% Reduce(cbind, .)
} else {
  metadata = metadatas %>% Reduce(rbind, .) 
  colnames(metadata)[ncol(metadata)] = 
    paste0("doubletProportion_", doubletProportion)
  metadata = metadata[, ncol(metadata), drop=FALSE]
  metadata[[grep("pANN_", colnames(metadata), value=TRUE)]]=NULL
}

print("Add doublet information to the metadata of the Seurat object ... ")
dat_filtered@meta.data = 
  cbind(dat_filtered@meta.data, metadata[rownames(dat_filtered@meta.data),])

print("Doublet QC ...")
doubletCols = grep("doubletProportion_", colnames(metadata), value=TRUE)
lapply(doubletCols, function(doubletCol) {
  pdf(paste0(OUTDIR, "/Doublets/coreQC_", doubletCol, ".pdf"), h=7.5, w=12)
    print(VlnPlot(dat_filtered, features=features_qc, 
          ncol=3, pt.size=0, split.by=doubletCol, group.by="Tissue"))
  dev.off()
  print("Doublet spatial QC ... ")
  Idents(dat_filtered) = dataset
  pdf(paste0(OUTDIR, "/Doublets/Doublets_spatialQC_", doubletCol, ".pdf"), h=10, w=10)
    print(ImageDimPlot(dat_filtered, fov = dataset, cols = "red", 
      cells = rownames(dat_filtered@meta.data)[dat_filtered@meta.data[[doubletCol]]=="Doublet"]) + 
        ggtitle(doubletCol))
  dev.off()
})


# dat_filtered = subset(dat_filtered, 
#   cells=rownames(dat_filtered@meta.data)[metadata$doublet_finder=="Singlet"])
print(dat_filtered)
print("Save the seurat object with doublet inforamtion...")
saveRDS(dat_filtered, paste0(OUTDIR, "/", dataset, "_filtered.RDS"))


print("QC log after remove doublets ...")
QClog = capture.output({
    cat("\n")
    print(dat_filtered@meta.data[, grep("doubletProportion_", 
      colnames(dat_filtered@meta.data), value=TRUE), drop=FALSE] %>%
          apply(., 2, table))
})

writeLines(QClog, paste0(OUTDIR, "/Doublets/QClog.txt"))


print("Save sessionInfo ... ")
out = capture.output(sessionInfo())
writeLines(out, paste0(OUTDIR, "/Doublets/sessionInfo.txt"))


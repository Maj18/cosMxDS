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
minSolidity = args_list[["minSolidity"]]
minArea.um2 = args_list[["minArea.um2"]]
minnCount_RNA = args_list[["minnCount_RNA"]]
INDIR = args_list[["INDIR"]]
doubletProportion = as.numeric(args_list[["doubletProportion"]])
print(minnCount_RNA)
# dataset = "Eleni_Female"
# outdir = "../results/QC/"
# minSolidity = 0.5
# minArea.um2 = 20
# minnCount_RNA = 25
# doubletProportion = 0.074
# INDIR = "/cfs/klemming/projects/supr/naiss2026-4-253/data/Shared_Folder/"




rmDoublets = function(data.filt, doubletProportion= 0.074) {
  # before doing doublet detection we need to run scaling, variable gene selection and PCA, as well as UMAP for visualization. 
  # These steps will be explored in more detail in coming exercises.
  DefaultAssay(data.filt) = "RNA"
  data.filt = NormalizeData(data.filt)
  data.filt = ScaleData(data.filt)
  data.filt = FindVariableFeatures(data.filt)
  data.filt = RunPCA(data.filt, verbose = F, npcs = 100)

  # Then we run doubletFinder, selecting first 10 PCs and a `pK` value of 0.9. To optimize the parameters, you can run the `paramSweep` function in the package.
  # set seed
  # set seed
  set.seed(8)

  # The expected doublet rate as suggested by 10X
  # The percentage of droplets containing two or more cells is influenced by the quantity of cells loaded into the sequencing machine.
  # ref: https://uofuhealth.utah.edu/huntsman/shared-resources/gcb/htg/single-cell/genomics-10x
  cells.nr = ncol(data.filt)
  nExp.pois = round(cells.nr * doubletProportion)
  # nExp.pois

  # pK identification (no ground-truth)
  # pK: neighborhood size, i.e. the number of neighbors (in percentage) to consider when calculating pANN
  sweep.list = invisible(paramSweep(data.filt, PCs = 1:10, num.cores = 1))
  sweep.stats = summarizeSweep(sweep.list)
  pdf("find.pK.plots.pdf")
  bcmvn = find.pK(sweep.stats)
  dev.off()

  # Optimal pK are selected by maximizing BCmvn (https://www.sciencedirect.com/science/article/pii/S2405471219300730)
  bcmvn.max = bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk = bcmvn.max$pK
  optimal.pk = as.numeric(levels(optimal.pk))[optimal.pk]
  # optimal.pk

  # run DoubletFinder
  # pN: the proportion of artificial doublet cells to generate (default: 0.25)
  data.filt = doubletFinder(seu = data.filt,
                          PCs = 1:10,
                          pK = optimal.pk,
                          nExp = nExp.pois)
  metadata = data.filt@meta.data
  colnames(metadata)[ncol(metadata)] = "doublet_finder"

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
dat = readRDS(paste0(INDIR, "/seuratObject_", dataset, ".RDS"))
dat

# Kmean clustering to define tissues
km = kmeans(dat@meta.data[,c("x_slide_mm", "y_slide_mm")], centers = 4, nstart = 25)
dat@meta.data$coord_cluster = km$cluster
temp = dat@meta.data %>% group_by(coord_cluster) %>% 
  summarize(meanx = mean(x_slide_mm), meany=mean(y_slide_mm)) %>%
  arrange(meanx, meany) %>% 
  mutate(Tissue=c(paste0("Injury", 1:2), paste0("Normal", 1:2)))
mapping = setNames(temp$Tissue, as.character(temp$coord_cluster))
pdf(paste0(OUTDIR, "/Tissue.pdf"), h=10, w=10)
  Idents(dat) = dat@meta.data$coord_cluster
  print(ImageDimPlot(dat, fov = dataset))
  dat = RenameIdents(dat,
                    `1` = mapping["1"],
                    `2` = mapping["2"],
                    `3` = mapping["3"],
                    `4` = mapping["4"])
  dat@meta.data$Tissue = Idents(dat)
  print(ImageDimPlot(dat, fov = dataset))
dev.off()
rm(mapping, temp)

# FOV summary 
# dat@meta.data = dat@meta.data %>% as.data.frame() %>%
#   group_by(fov) %>%
#   summarize(nCount_RNA_FOVmean = mean(nCount_RNA)) %>%
#   right_join(dat@meta.data)


# Before filtering
plotQC = function(dat, outDIR=paste0(OUTDIR, "/beforeFiltering/")) {
  ## core QC
  features_qc = c(
    "nCount_RNA",
    "nFeature_RNA",
    "Area.um2",
    "Circularity", 
    "Solidity",
    "Mean.rRNA",
    "complexity",
    "propNegative",
    "percOfDataFromError"
  )
  dat@meta.data$Tissue = 
    factor(dat@meta.data$Tissue, levels = c(paste0("Normal", 1:2), paste0("Injury", 1:2)))
  print("Violin plot...")
  pdf(paste0(outDIR, "/coreQC.pdf"), h=7.5, w=12)
  print(VlnPlot(dat, features = features_qc, ncol = 3, pt.size=0, group.by="Tissue"))
  dev.off()

  ## Scatter relationships
  print("Scatter plot ...")
  Idents(dat) = dataset
  pdf(paste0(outDIR, "/QC_scatterplot.pdf"), h=3, w=7)
  p1 = FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "propNegative")
  p2 = FeatureScatter(dat, feature1 = "Area.um2", feature2 = "nCount_RNA")
  print(p1+p2)
  dev.off()

  ## Spatial QC
  print("Spatial QC ...")
  dir.create(paste0(outDIR, "/spatialQC"), recursive=T, showWarnings=FALSE)
  ### scale some of the metadata
  meta.data = dat@meta.data
  for (f in features_qc[c(1:6)]) {
    meta.data[[paste0("log10", f)]] = meta.data[[f]]
    meta.data[[paste0("log10", f)]] = log10(meta.data[[f]]+1)
  }
  
  dat@meta.data = meta.data
  lapply(c(paste0("log10",features_qc[c(1:6)]), features_qc[7:9], "qcCellsPassed"), function(f) {
    pdf(paste0(outDIR, "/spatialQC/coreQC_spatial_", f, ".pdf"), h=5, w=5)
    print(ImageFeaturePlot(
      dat, features = f, size=0.1, border.color="white", border.size=0.1))
    dev.off()
  })

  return(dat)
}

outDIR=paste0(OUTDIR, "/beforeFiltering/")
dir.create(outDIR, recursive=T, showWarnings=FALSE)
dat = plotQC(dat, outDIR=outDIR)
## Extra QC plots
print("extra QC plots ...")
pdf(paste0(outDIR, "/extraQC.pdf"), h=10, w=10)
  print(ImageDimPlot(dat, fov = dataset, cols = "red", 
    cells = WhichCells(dat, expression=Area.um2 < minArea.um2))+ggtitle(paste0("Area.um2<",minArea.um2)))
  print(ImageDimPlot(dat, fov = dataset, cols = "red", 
    cells = WhichCells(dat, expression=Solidity < minSolidity))+ggtitle(paste0("Solidity<",minSolidity)))
  print(ImageDimPlot(dat, fov = dataset, cols = "red", 
    cells = WhichCells(dat, expression=nCount_RNA < minnCount_RNA))+ggtitle(paste0("nCount_RNA<",minnCount_RNA)))
dev.off()



# QC log
print("QC log before filtering ...")
QClogb = capture.output({
  cat("Before filtering...\n")
  features_qc = c(
    "nCount_RNA",
    "nFeature_RNA",
    "Area.um2",
    "Circularity", 
    "Solidity",
    "propNegative",
    "complexity",
    "percOfDataFromError",
    "Mean.rRNA"
  )
  dat@meta.data[, features_qc] %>% apply(., 2, summary)
  flags = c("qcFlagsCellCounts",
  "qcFlagsCellPropNeg",
  "qcFlagsCellComplex",
  "qcFlagsCellArea",
  "qcCellsFlagged",
  "qcCellsPassed",
  "qcFlagsFOV")
  grep("qc", colnames(dat@meta.data), value=T)
  print(dat@meta.data[, flags] %>%
    as.data.frame() %>%
    mutate(qcCellsFlagged = ifelse(qcCellsFlagged, "Fail", "Pass")) %>%
    mutate(qcCellsPassed = ifelse(qcCellsPassed, "Pass", "Fail")) %>% 
    apply(., 2, table))
  cat("\n")
  dat@meta.data[, features_qc] %>% apply(., 2, summary)
})

print("QC log filtering ...")
QClogf = capture.output({
  cat ("Keep a copy of failed cells:\n")
  rownames(dat@meta.data)[
    !(dat@meta.data$qcFlagsCellCounts == "Pass" &
    dat@meta.data$qcFlagsCellPropNeg == "Pass" &
    dat@meta.data$qcFlagsCellComplex == "Pass" &
    dat@meta.data$qcFlagsCellArea == "Pass" &
    !dat@meta.data$qcCellsFlagged &
    dat@meta.data$qcCellsPassed &
    dat@meta.data$qcFlagsFOV == "Pass" &
    dat@meta.data$Solidity > minSolidity &
    dat@meta.data$Area.um2 > minArea.um2 &
    dat@meta.data$nCount_RNA > minnCount_RNA)
  ]
})

print("Filter cells")
passedCells = rownames(dat@meta.data)[
  (dat@meta.data$qcFlagsCellCounts == "Pass" & # Cell QC
  dat@meta.data$qcFlagsCellPropNeg == "Pass" & # Cell QC
  dat@meta.data$qcFlagsCellComplex == "Pass" & # Cell QC
  dat@meta.data$qcFlagsCellArea == "Pass" & # Cell segmentation QC
  !dat@meta.data$qcCellsFlagged & # Cell QC
  dat@meta.data$qcCellsPassed &  # Cell QC
  dat@meta.data$qcFlagsFOV == "Pass" & # FOV QC
  dat@meta.data$Solidity > minSolidity & # Cell segmentation QC
  dat@meta.data$Area.um2 > minArea.um2 & # Cell segmentation QC
  dat@meta.data$nCount_RNA > minnCount_RNA) # Cell QC
]

print("Only keep genes that non-0 counts in more than 5 cells ...")
passedGenes = 
  rownames(dat@assays$RNA@counts)[rowSums(dat@assays$RNA@counts>0) > 5]
dat = UpdateSeuratObject(dat)
dat_filtered = subset(dat, cells=passedCells, features=passedGenes)
# saveRDS(dat_filtered, paste0(OUTDIR, "/", dataset, "_filtered.RDS"))
# dat_filtered = readRDS(paste0(OUTDIR, "/", dataset, "_filtered.RDS"))
# dat_filtered$smallCells = ifelse(dat_filtered@meta.data$Area.um2 < 20, "small", "large")

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
metadata = rmDoublets(data.filt=dat_filtered, doubletProportion= doubletProportion)
pdf(paste0(outDIR, "/Doublets_QC.pdf"), h=10, w=10)
print(ImageDimPlot(dat_filtered, fov = dataset, cols = "red", 
      cells = rownames(dat_filtered@meta.data)[metadata$doublet_finder=="Doublet"]) + 
        ggtitle("Doublets"))
dev.off()
dat_filtered@meta.data = metadata
print("Remove doublets ...")
dat_filtered = subset(dat_filtered, 
  cells=rownames(dat_filtered@meta.data)[metadata$doublet_finder=="Singlet"])
saveRDS(dat_filtered, paste0(OUTDIR, "/", dataset, "_filtered.RDS"))



print("QC log after filtering ...")
QCloga = capture.output({
  cat("After filtering...\n")
  dat_filtered@meta.data[, flags] %>%
    as.data.frame() %>%
    mutate(qcCellsFlagged = ifelse(qcCellsFlagged, "Fail", "Pass")) %>%
    mutate(qcCellsPassed = ifelse(qcCellsPassed, "Pass", "Fail")) %>% 
    apply(., 2, table)
  cat("\n")
  dat_filtered@meta.data[, features_qc] %>% apply(., 2, summary)
})

QClog = c(QClogb, "\n\n", QCloga, "\n\n", QClogf)
writeLines(QClog, paste0(OUTDIR, "/QClog.txt"))



# Plot QC after filteirng:
print("Plot QC after filteirng ...")
outDIR2=paste0(OUTDIR, "/afterFiltering/")
dir.create(outDIR2, recursive=T, showWarnings=FALSE)
plotQC(dat_filtered, outDIR=outDIR2)
quantile(dat_filtered@meta.data$nCount_RNA, probs=0.005)

out = capture.output(sessionInfo())
writeLines(out, paste0(OUTDIR, "/sessionInfo.txt"))



#!/opt/conda/bin/R
args = commandArgs(trailingOnly = TRUE)
# Parse key=value pairs
args_list = list()
for (arg in args) {
  kv = strsplit(arg, "=")[[1]]
  args_list[[kv[1]]] = kv[2]
}
dataset = args_list[["batch"]]
outdir = args_list[["outdir"]]
minSolidity = args_list[["minSolidity"]]
minArea.um2 = args_list[["minArea.um2"]]
minnCount_RNA = args_list[["minnCount_RNA"]]
# dataset = "Eleni_Female"
# outdir = "../results/QC/"
# minSolidity = 0.5
# minArea.um2 = 20
# minnCount_RNA = 25

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


OUTDIR = paste0(outdir, "/",dataset, "/")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dat = readRDS(paste0("../data/Shared_Folder/seuratObject_", dataset, ".RDS"))
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
  pdf(paste0(outDIR, "/coreQC.pdf"), h=7.5, w=12)
  print(VlnPlot(dat, features = features_qc, ncol = 3, pt.size=0, group.by="Tissue"))
  dev.off()

  ## Scatter relationships
  pdf(paste0(outDIR, "/QC_scatterplot.pdf"), h=3, w=7)
  p1 = FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "propNegative")
  p2 = FeatureScatter(dat, feature1 = "Area.um2", feature2 = "nCount_RNA")
  print(p1+p2)
  dev.off()

  ## Spatial QC
  dir.create(paste0(outDIR, "/spatialQC"), recursive=T, showWarnings=FALSE)
  ### scale some of the metadata
  meta.data = dat@meta.data
  for (f in features_qc[c(1:6)]) {
    meta.data[[paste0("log10", f)]] = meta.data[[f]]
    meta.data[[paste0("log10", f)]] = log10(meta.data[[f]]+1)
  }
  
  dat@meta.data = meta.data
  lapply(c(paste0("log10",features_qc[1:6]), features_qc[7:9], "qcCellsPassed"), function(f) {
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
pdf(paste0(outDIR, "/extraQC.pdf"), h=10, w=10)
  print(ImageDimPlot(dat, fov = dataset, cols = "red", 
    cells = WhichCells(dat_filtered, expression=Area.um2 < minArea.um2))+ggtitle(paste0("Area.um2<",minArea.um2)))
  print(ImageDimPlot(dat, fov = dataset, cols = "red", 
    cells = WhichCells(dat_filtered, expression=Solidity < minSolidity))+ggtitle(paste0("Solidity<",minSolidity)))
  print(ImageDimPlot(dat, fov = dataset, cols = "red", 
    cells = WhichCells(dat_filtered, expression=nCount_RNA < minnCount_RNA))+ggtitle(paste0("nCount_RNA<",minnCount_RNA)))
dev.off()



# QC log
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
    Solidity > minSolidity &
    Area.um2 > minArea.um2 &
    nCount_RNA > minnCount_RN)
  ]
})

passedCells = rownames(dat@meta.data)[
  (dat@meta.data$qcFlagsCellCounts == "Pass" &
  dat@meta.data$qcFlagsCellPropNeg == "Pass" &
  dat@meta.data$qcFlagsCellComplex == "Pass" &
  dat@meta.data$qcFlagsCellArea == "Pass" &
  !dat@meta.data$qcCellsFlagged &
  dat@meta.data$qcCellsPassed &
  dat@meta.data$qcFlagsFOV == "Pass" &
  Solidity > minSolidity &
  Area.um2 > minArea.um2 &
  nCount_RNA > minnCount_RN)
]
dat = UpdateSeuratObject(dat)
dat_filtered = subset(dat, cells=passedCells)
saveRDS(dat_filtered, paste0(OUTDIR, "/", dataset, "_filtered.RDS"))
# dat_filtered = readRDS(paste0(OUTDIR, "/", dataset, "_filtered.RDS"))
# dat_filtered$smallCells = ifelse(dat_filtered@meta.data$Area.um2 < 20, "small", "large")

QCloga = capture.output({
  cat("After filtering...\n")
  dat_filtered@meta.data[, flags] %>%
    as.data.frame() %>%
    mutate(qcCellsFlagged = ifelse(qcCellsFlagged, "Fail", "Pass")) %>%
    mutate(qcCellsPassed = ifelse(qcCellsPassed, "Pass", "Fail")) %>% 
    apply(., 2, table)
  dat_filtered@meta.data[, features_qc] %>% apply(., 2, summary)
})

QClog = c(QClogb, "\n\n", QCloga, "\n\n", QClogf)
writeLines(QClog, paste0(OUTDIR, "/QClog.txt"))



# Plot QC after filteirng:
outDIR2=paste0(OUTDIR, "/afterFiltering/")
dir.create(outDIR2, recursive=T, showWarnings=FALSE)
plotQC(dat, outDIR=outDIR2)
quantile(dat@meta.data$nCount_RNA, probs=0.005)

out = capture.output(sessionInfo())
writeLines(out, paste0(OUTDIR, "/sessionInfo.txt"))

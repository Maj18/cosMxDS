
#!/opt/conda/bin/R
args = commandArgs(trailingOnly = TRUE)
# Parse key=value pairs
args_list = list()
for (arg in args) {
  kv = strsplit(arg, "=")[[1]]
  args_list[[kv[1]]] = kv[2]
}
batch = args_list[["batch"]]
outdir = args_list[["outdir"]]
# dataset = "Eleni_Female"
# outdir = "../results/QC/"

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


OUTDIR = paste0(outdir, "/",dataset, "/")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dat = readRDS(paste0("../data/Shared_Folder/seuratObject_", dataset, ".RDS"))
dat

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
  pdf(paste0(outDIR, "/coreQC.pdf"), h=7.5, w=7.5)
  print(VlnPlot(dat, features = features_qc, ncol = 3, pt.size=0))
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
    dat@meta.data$qcFlagsFOV == "Pass")
  ]
})

passedCells = rownames(dat@meta.data)[
  (dat@meta.data$qcFlagsCellCounts == "Pass" &
  dat@meta.data$qcFlagsCellPropNeg == "Pass" &
  dat@meta.data$qcFlagsCellComplex == "Pass" &
  dat@meta.data$qcFlagsCellArea == "Pass" &
  !dat@meta.data$qcCellsFlagged &
  dat@meta.data$qcCellsPassed &
  dat@meta.data$qcFlagsFOV == "Pass")
]
dat = UpdateSeuratObject(dat)
dat_filtered = subset(dat, cells=passedCells)
saveRDS(dat_filtered, paste0(OUTDIR, "/", dataset, "_filtered.RDS"))
# dat_filtered = readRDS(paste0(OUTDIR, "/", dataset, "_filtered.RDS"))

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

out = capture.output(sessionInfo())
writeLines(out, paste0(OUTDIR, "/sessionInfo.txt"))


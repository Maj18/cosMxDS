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

batches = args_list[["batches"]] %>% strsplit(.,",") %>% .[[1]]
INFILES = args_list[["INFILES"]] %>% strsplit(.,",") %>% .[[1]]
outdir = paste0(args_list[["outdir"]], "/QC/")
outdir %>% dir.create(., recursive=TRUE, showWarnings=FALSE)

print(batches)
print(INFILES)
print(outdir)


# batches = "Eleni_Female,Eleni_Male"
# INFILES = "../results/QC/Eleni_Female/Eleni_Female_filtered.RDS,../results/QC/Eleni_Female/Eleni_Female_filtered.RDS"
# outdir = "../results/QC/"


print("Load packages...")
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

print("Import seurat objects from all batches ...")
dat.list = lapply(INFILES, function(INFILE) {
  readRDS(INFILE)
}) 

print("Merge two batches")
dat = merge(dat.list[[1]], dat.list[[2]], add.cell.ids = gsub("Eleni_", "", batches))
rm(dat.list)
gc()

print(summary cells by batch and tissue ...)
plotQC = function(dat, outDIR=outdir) {
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
    pdf(paste0(outDIR, "/coreQC_all.pdf"), h=7.5, w=12)
    print(VlnPlot(dat, features = features_qc, ncol = 3, pt.size=0, group.by="Tissue"))
    dev.off()
    
    summary = dat.list[[1]]@meta.data %>% as.data.frame() %>% 
      group_by(Run_Tissue_name, Tissue) %>% summarize(
        CellNr = n(),
        mediannCount_RNA = median(nCount_RNA),
        mediannFeature_RNA = median(nFeature_RNA),
        medianArea.um2 = median(Area.um2),
        medianCircularity = median(Circularity),
        medianSolidity = median(Solidity),
        medianMean.rRNA = median(Mean.rRNA),
        mediancomplexity = median(complexity),
        medianpropNegative = median(propNegative),
        medianpercOfDataFromError = median(percOfDataFromError))
    write.table(summary, paste0(outdir, "Summary_all.csv"), quote=FALSE, row.names=FALSE)
}

plotQC(dat, outDIR=outdir)

saveRDS(dat, function(outdir, "All_filtered.RDS"))

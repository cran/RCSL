## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE, include = FALSE----
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png',
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("QinglinMei/RCSL")

## ---- results="hide"----------------------------------------------------------
library(RCSL)
library(SingleCellExperiment)
library(ggplot2)
library(igraph)
library(umap)

## -----------------------------------------------------------------------------
head(ann)
yan[1:3, 1:3]
origData <- yan
label <- ann$cell_type1

## ---- cache=TRUE--------------------------------------------------------------
data <- log2(as.matrix(origData) + 1)
gfData <- GenesFilter(data)

## ---- cache=TRUE--------------------------------------------------------------
resSimS <- SimS(gfData)

## ---- cache=TRUE--------------------------------------------------------------
Estimated_C <- EstClusters(resSimS$drData,resSimS$S)

## ---- cache=TRUE--------------------------------------------------------------
resBDSM <- BDSM(resSimS$S, Estimated_C)

## ---- cache=TRUE--------------------------------------------------------------
ARI_RCSL <- igraph::compare(resBDSM$y, label, method = "adjusted.rand")

## ---- cache=TRUE--------------------------------------------------------------
DataName <- "Yan"
res_TrajecAnalysis <- TrajectoryAnalysis(gfData, resSimS$drData, resSimS$S,
                                         clustRes = resBDSM$y, TrueLabel = label, 
                                         startPoint = 1, dataName = DataName)

## ---- cache=TRUE--------------------------------------------------------------
res_TrajecAnalysis$MSTPlot

## ---- cache=TRUE--------------------------------------------------------------
res_TrajecAnalysis$PseudoTimePlot

## ---- cache=TRUE--------------------------------------------------------------
res_TrajecAnalysis$TrajectoryPlot


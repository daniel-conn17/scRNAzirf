#' Run SCENIC using zirf VIMs
#'
#' @param scenicOptions Options to set up SCENIC runs
#' @param weightMatrix_list List of VIMs from ZIRFs, designed to take in the output
#'        of zirf_genie3.
#' @param exprMatrix Normalized expression matrix.
#'
#' @param log_scale Should log of exprMatrix be used?
#' @return As in SCENIC, the output is written in the folders 'int' and 'output'
#' @export
run_SCENIC <- function(scenicOptions, weightMatrix_list, exprMatrix, log_scale = FALSE){
  linkList_list <- list()
  weightMatrix_list <- lapply(weightMatrix_list, as.matrix)
  for(i in 1:length(weightMatrix_list)) {
    # If interrupted or run by pieces:
    # fileName <- gsub(".Rds$", paste0("_part_", i,".Rds"), getIntName(scenicOptions, "genie3wm"))
    # weightMatrix <- readRDS(fileName)
    weightMatrix <- weightMatrix_list[[i]]
    linkList_list[[i]] <- GENIE3::getLinkList(weightMatrix, threshold=SCENIC::getSettings(scenicOptions, "modules/weightThreshold"))
    linkList_list[[i]]$regulatoryGene <- as.character(linkList_list[[i]]$regulatoryGene)
    linkList_list[[i]]$targetGene <- as.character(linkList_list[[i]]$targetGene)
    linkList_list[[i]]$weight <- as.numeric(linkList_list[[i]]$weight)
  }

  #rm(weightMatrices)

  # Merge
  linkList <- do.call("rbind", linkList_list)
  colnames(linkList) <- c("TF", "Target", "weight")

  #linkList[, "TF"] <- fix_Xgenes(linkList[, "TF"])
  #linkList[, "Target"] <- fix_Xgenes(linkList[, "Target"])

  # Order by weight
  linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
  saveRDS(linkList, file=SCENIC::getIntName(scenicOptions, "genie3ll"))

  #Next run the following code from the following link:
  #https://github.com/aertslab/SCENIC/blob/master/vignettes/detailedStep_1_coexNetwork2modules.Rmd
  linkList <- SCENIC::loadInt(scenicOptions, "genie3ll")
  if(!all(colnames(linkList) == c("TF", "Target", "weight"))) stop('The link list colnames should be "TF", "Target", "weight"')

  table(linkList$Target)
  table(linkList$Target)
  uniquePairs <- nrow(unique(linkList[,c("TF", "Target")]))
  if(uniquePairs < nrow(linkList))
    stop("There are duplicated regulator-target (gene id/name) pairs in the input link list.")

  msg <- paste0(format(Sys.time(), "%H:%M"), "\tCreating TF modules")
  if(SCENIC::getSettings(scenicOptions, "verbose")) message(msg)
  print(stats::quantile(linkList$weight, probs=c(0.75, 0.90)))
  .openDev <- function(fileName, devType, ...) {
    if(devType=="pdf")
      grDevices::pdf(paste0(fileName, ".pdf"), ...)

    if(devType=="png")
      grDevices::png(paste0(fileName, ".png"), ...)

    if(devType=="cairo_pfd") # similar to Cairo::CairoPDF?
      grDevices::cairo_pdf(paste0(fileName, ".pdf"), ...)
  }
  .openDev(fileName=SCENIC::getIntName(scenicOptions, "genie3weighPlot"),
           devType=SCENIC::getSettings(scenicOptions, "devType"))
  plot(linkList$weight[1:1000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
       ylab="Weight", xlab="Links sorted decreasingly")
  graphics::abline(h=0.001, col="blue") # Threshold
  #sum(linkList$weight>0.001)/nrow(linkList)


  SCENIC::runCorrelation(exprMatrix, scenicOptions)
  if(!log_scale) {
    exprMatrix_log <- log2(exprMatrix+1)
  }


  cntPairs <- table(table(linkList[,"TF"],linkList[,"Target"]))
  if(any(names(cntPairs)>1))
    stop("There are duplicated regulator-target (gene id/name) pairs in the input link list.")

  SCENIC::runSCENIC_1_coexNetwork2modules(scenicOptions)
  SCENIC::runSCENIC_2_createRegulons(scenicOptions,
                             coexMethod = c("top50", "top5perTarget", "top10perTarget", "top50perTarget"))
  SCENIC::runSCENIC_3_scoreCells(scenicOptions, exprMatrix)
  SCENIC::runSCENIC_4_aucell_binarize(scenicOptions, exprMatrix)
}

### Plot impulse fits

#' @include srcImpulseDE2_evalImpulse.R
#' @include srcImpulseDE2_classImpulseDE2Object.R
NULL

#' Plots the impulse fits and data
#' 
#' Plots the impulse fits and data to pdf and return a list of ggplots. Points
#' are size-factor and covariate- adjusted data. Consider using boolSimplePlot=TRUE
#' if the plot seems to crowded.
#' 
#' @seealso Called by \link{computeModelFits}.
#' 
#' @param vecGeneIDs (string vector) [Default NULL]
#' Gene names to be plotted. Must be in rownames of 
#' objectImpulseDE2@matCountDataProc.
#' Supply either vecGeneIDs or scaNTopIDs.
#' @param scaNTopIDs (int) [Default NULL]
#' Number of top differentially expressed (by q-value) genes to 
#' be plotted
#' Supply either vecGeneIDs or scaNTopIDs.
#' @param objectImpulseDE2 (ImpulseDE2 object)
#' Object previously fitted to be used for plotting.
#' @param boolCaseCtrl (bool) Whether to create case-ctrl plot.
#' @param dirOut (dir) [Default NULL]
#' Directory into which pdf is printed.
#' @param strFileName (str) [Default 'ImpulseDE2_Trajectories.pdf']
#' File name of pdf with plots.
#' @param boolMultiplePlotsPerPage (bool) [Default TRUE]
#' Whether to create grid with multiple plots on each page of pdf.
#' @param boolSimplePlot (bool) [Default TRUE]
#' Whether to omit batch structure in plotting of model fits
#' and only plot fit to first batch/all data (if no confounders were given).
#' This strongly simplifies plots and is recommended e.g. for case-ctrl data.
#' @param vecRefPval (vector length vecGeneIDs) [Default NULL]
#' P/Q-values to be displayed alongside ImpulseDE2 q-value
#' for differential expression in plot titles.
#' @param strNameRefMethod (str) [Default NULL]
#' Name of reference method used to generate vecRefPval.
#' Mentioned in plot titles.
#' 
#' @return lsgplotsID (gplot list length vecGeneIDs)
#' List of gplots for IDs in vecGeneIDs. This is secondary 
#' output next to the .pdf and can be used to extract 
#' single plots or assemble plots differently.
#' 
#' @examples
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 0,
#' scaNImp          = 40,
#' scaNLin          = 20,
#' scaNSig          = 40)
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' vecConfounders  = NULL,
#' boolIdentifyTransients = FALSE,
#' scaNProc        = 1 )
#' lsgplotsID <- plotGenes(
#' scaNTopIDs=5,
#' objectImpulseDE2=objectImpulseDE2,
#' boolCaseCtrl=FALSE,
#' boolMultiplePlotsPerPage=TRUE,
#' boolSimplePlot=FALSE)
#' lsgplotsID[[1]]
#' 
#' @author David Sebastian Fischer
#' 
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' 
#' @export
plotGenes <- function(
    vecGeneIDs = NULL, scaNTopIDs = NULL, objectImpulseDE2, 
    boolCaseCtrl, dirOut = NULL, strFileName = "ImpulseDE2_Trajectories.pdf", 
    boolMultiplePlotsPerPage = TRUE, boolSimplePlot = FALSE, vecRefPval = NULL, 
    strNameRefMethod = NULL) {
    
    dfAnnot <- get_dfDESeqAnnotationProc(obj=objectImpulseDE2)
    # Set graphical parameters
    scaNPlotsPerPage <- 4
    # Colour-blind palette
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7")
    
    # Check input
    if (is.null(get_lsModelFits(objectImpulseDE2))) {
        error(paste0("objectImpulseDE2 does not contain model fits. ",
                     "Run ImpulseDE2_main first."))
    }
    if (is.null(vecGeneIDs) & is.null(scaNTopIDs)) {
        stop("Supply either vecGeneIDs or scaNTopIDs.")
    }
    if (!is.null(vecGeneIDs) & !is.null(scaNTopIDs)) { 
        stop("Only one of the two: vecGeneIDs or scaNTopIDs.")
    }
    if (is.null(get_vecCovFactor(objectImpulseDE2))) {
        print("Setting boolSimplePlot=TRUE as no categorical covariates were found.")
        boolSimplePlot <- TRUE
    }
    if (!is.null(dirOut) && !file.exists(dirOut)) {
        stop("Output directory dirOut not available.")
    }
    if (!is.null(vecRefPval) && (names(vecRefPval) != vecGeneIDs)) {
        stop("Names of vecRefPval have to be IDs from vecGeneIDs.")
    }
    if (!boolSimplePlot & boolCaseCtrl & boolMultiplePlotsPerPage) {
        warning(paste0("Plots are likely overloaded. ",
                       "Consider switching to boolSimplePlot=TRUE ",
                       "or boolMultiplePlotsPerPage=FALSE."))
    }
    # Use top scaNTopIDs IDs if no specific IDs were supplied.
    if (is.null(vecGeneIDs)) 
        vecGeneIDs <- objectImpulseDE2$dfImpulseDE2Results[
            with(objectImpulseDE2$dfImpulseDE2Results, 
                 order(padj)), ]$Gene[1:scaNTopIDs]
    scaNIDs <- length(vecGeneIDs)
    lsgplotsID <- list()

    # Get size-factor, covariate-adjusted, sample-level values for points 
    # Gene x sample data frame 
    dfGeneSample <- computeModelFits(objectImpulseDE2)$sample

    for (id in vecGeneIDs) {
        # Sequence of times to fit 
        vecTimePointsFit <- seq(
            min(get_dfDESeqAnnotationProc(obj=objectImpulseDE2)$Time), 
                                max(get_dfDESeqAnnotationProc(obj=objectImpulseDE2)$Time), 
                                length.out = 500)

        # Values from fitted impulse model - "case"
        vecCaseImpulseParam <- get_lsModelFits(
            obj=objectImpulseDE2)$case[[id]]$lsImpulseFit$vecImpulseParam
        vecCaseImpulseValue <- evalImpulse_comp(
            vecImpulseParam = vecCaseImpulseParam, 
            vecTimepoints = vecTimePointsFit)

        # Get size-factor, covariate-adjusted, sample-level values for points 
        dfPoint = data.frame(t(dfGeneSample[id,]))
        rownames(dfPoint) = colnames(dfGeneSample)
        colnames(dfPoint) = "y"
        dfPoint$time = dfAnnot$Time[match(rownames(dfPoint), dfAnnot$Sample)]

        title = paste0(id, ": log10 padj ", 
                    round(log10(objectImpulseDE2$dfImpulseDE2Results[id, ]$padj)))

        # Get batch adjustments, if applicable 
        if(is.null(get_vecCovFactor(objectImpulseDE2)) | boolSimplePlot){
            # If there is more than one, pick the one with the fewest levels 
            vecCovFactor = get_vecCovFactor(objectImpulseDE2)
            if(length(vecCovFactor) > 1){
                vecNumLvls = as.numeric(unlist(lapply(vecCovFactor, function(x){
                    length(unique(dfAnnot[,x])) - 1
                    })))
                names(vecNumLvls) = vecCovFactor
                covFactor = names(vecNumLvls)[which.min(vecNumLvls)]
            }else{
                covFactor = vecCovFactor
            }
            # Get scaling factors 
            vecFactorCoef = get_lsModelFits(obj=obj)$case[[x]]$lsImpulseFit$vecCorrectionFactors
            vecFactorCoef = vecFactorCoef[grepl(covFactor, names(vecFactorCoef))]
            lsvecFactorCoef = list(case=vecFactorCoef)
            if(boolCaseCtrl){
                # Repeat it for "control" and "combined"
                for(c in c("control","combined")){
                    vecFC = get_lsModelFits(obj=obj)[[c]][[x]]$lsImpulseFit$vecCorrectionFactors
                    vecFC = vecFC[grepl(covFactor, names(vecFC))]
                    lsvecFactorCoef[[c]] = vecFC
                }
            }else{
                lsvecFactorCoef$control = NULL
                lsvecFactorCoef$combined = NULL
            }
        }else{
            lsvecFactorCoef = list(case = NULL, control = NULL, combined = NULL)
        }

        if(boolCaseCtrl){
            # Three solid lines - one for "case", one for "control", one for "combined"
            # Values from fitted impulse model - "control"
            vecCtrlImpulseParam <- get_lsModelFits(
                obj=objectImpulseDE2)$control[[id]]$lsImpulseFit$vecImpulseParam
            vecCtrlImpulseValue <- evalImpulse_comp(
                vecImpulseParam = vecCtrlImpulseParam, 
                vecTimepoints = vecTimePointsFit)
            # Values from fitted impulse model - "combined"
            vecCmbImpulseParam <- get_lsModelFits(
                obj=objectImpulseDE2)$combined[[id]]$lsImpulseFit$vecImpulseParam
            vecCmbImpulseValue <- evalImpulse_comp(
                vecImpulseParam = vecCmbImpulseParam, 
                vecTimepoints = vecTimePointsFit)
            dfCase = data.frame(time=vecTimePointsFit, y=vecCaseImpulseValue, condition="case")
            dfCtrl = data.frame(time=vecTimePointsFit, y=vecCtrlImpulseValue, condition="control")
            dfCmb = data.frame(time=vecTimePointsFit, y=vecCmbImpulseValue, condition="combined")
            dfTime = rbind(dfCase, dfCtrl, dfCmb)

            # Add "condition" to dfPoint (no "combined")
            dfPoint$condition = dfAnnot$Condition[match(rownames(dfPoint), dfAnnot$Sample)]

            if(is.null(lsvecFactorCoef$case)){

                # Start the plot with points and solid lines; colour = condition, not batch 
                plotGene = ggplot() +
                    geom_point(data=dfPoint, aes(x=time, y=y, colour=condition)) +
                    geom_line(data=dfLine, aes(x=time, y=y, colour=condition)) +
                    theme_classic() + 
                    scale_colour_manual(values=c(case="blue", control="green3", combined="black"), name="Condition") +
                    labs(title=title, x='Time', y='Adjusted counts')
            }else{
                # Add non-solid lines for non-reference batches
                lsdfLine = list()
                i = 1
                for(c in c("case", "control", "combined")){
                    for(lvl in names(lsvecFactorCoef$c)){
                        dfTmp = dfLine[dfLine$condition == c]
                        dfTmp$y = dfTmp$y * lsvecFactorCoef[[lvl]]
                        dfTmp$batch = lvl
                        dfTmp$condition = c
                        lsdfLineCase[[i]] = dfTmp
                        i = i+1
                    }
                }
                dfLineBatch = do.call("rbind", lsdfLine) 

                # define reference "dfLine" level 
                all_levels = unique(paste0(covFactor, dfAnnot$covFactor))
                ref = all_levels[!all_levels %in% names(lsvecFactorCoef$combined)]
                dfLine$condition = ref
                ref_linetype = "solid"
                other_linetype = c("blank","dashed","dotted","dotdash","longdash","twodash")[1:min(vecNumLvls)]
                plotLines = c(ref_linetype, other_linetype)
                names(plotLines) = c(ref, unique(dfLineBatch$batch))

                plotGene = ggplot() +
                    geom_point(data=dfPoint, aes(x=time, y=y, colour=condition)) +
                    geom_line(data=dfLine, aes(x=time, y=y, colour=condition)) +
                    theme_classic() + 
                    scale_colour_manual(values=c(case="blue", control="green3", combined="black"), name="Condition") +
                    scale_linetype_manual(values=plotLines) +
                    labs(title=title, x='Time', y='Adjusted counts')
            }

        }else{
            if(is.null(lsvecFactorCoef$case)){
                # One solid line ("case")
                dfLine = data.frame(time=vecTimePointsFit, y=vecCaseImpulseValue)
                dfPoint = data.frame(t(dfGeneSample[id,]))
                rownames(dfPoint) = colnames(dfGeneSample)
                colnames(dfPoint) = "y"
                dfPoint$time = dfAnnot$Time[match(rownames(dfPoint), dfAnnot$Sample)]

                # Start the plot with points and solid lines
                plotGene = ggplot() +
                    geom_point(data=dfPoint, aes(x=time, y=y)) +
                    geom_line(data=dfLine, aes(x=time, y=y), colour='black') +
                    theme_classic() + 
                    labs(title=title, x='Time', y='Adjusted counts')
            }else{
                # Colour points by batch 
                # One solid line ("case")
                dfLine = data.frame(time=vecTimePointsFit, y=vecCaseImpulseValue)
                dfPoint = data.frame(t(dfGeneSample[id,]))
                rownames(dfPoint) = colnames(dfGeneSample)
                colnames(dfPoint) = "y"
                dfPoint$time = dfAnnot$Time[match(rownames(dfPoint), dfAnnot$Sample)]
                # Add batch to dfPoint
                dfPoint$batch = dfAnnot[,covFactor][match(rownames(dfPoint), dfAnnot$Sample)]
                dfPoint$batch = paste0(covFactor, dfPoint$batch)

                # Add to dfLine
                lsdfLine = list()
                i = 1
                for(lvl in names(lsvecFactorCoef$case)){
                    dfTmp = dfLine[dfLine$condition == "case"]
                    dfTmp$y = dfTmp$y * lsvecFactorCoef[[lvl]]
                    dfTmp$batch = lvl
                    lsdfLineCase[[i]] = dfTmp
                    i = i+1
                }
                dfLineBatch = do.call("rbind", lsdfLine) 
                
                # make the reference line black
                ref = unique(dfPoint$batch[!dfPoint$batch %in% dfLineBatch$batch])
                ref_color = "black"
                other_colors = cbPalette[1:min(vecNumLvls)]
                plotCols = c(ref_color, other_colors)
                names(plotCols) = c(ref, unique(dfLineBatch$batch))

                dfLine$batch = ref
                dfLine = data.frame(rbind(dfLine, dfLineBatch))

                # Add covFactor to dfPoint 
                # Start the plot with points and solid lines
                plotGene = ggplot() +
                    geom_point(data=dfPoint, aes(x=time, y=y, colour=batch)) +
                    geom_line(data=dfLine, aes(x=time, y=y, colour=batch)) +
                    theme_classic() + 
                    labs(title=title, x='Time', y='Adjusted counts') +
                    geom_colour_manual(values=plotCols) 
            }
        }
        lsgplotsID[[length(lsgplotsID) + 1]] <- gplotID
    }
    
    # Print to file
    if (!is.null(dirOut)) {
        dirFileOut <- paste0(dirOut, strFileName)
        print(paste0("Creating ", dirFileOut))
        graphics.off()
        if (boolMultiplePlotsPerPage) {
            pdf(dirFileOut)
            scaNPages <- scaNIDs%/%scaNPlotsPerPage
            if (scaNIDs%%scaNPlotsPerPage == 0) 
                scaNPages <- scaNPages - 1
            for (p in seq(0, scaNPages)) {
                if (p < scaNIDs%/%scaNPlotsPerPage) {
                    vecidxPlots <- seq((p * scaNPlotsPerPage + 1), 
                                       ((p + 1) * (scaNPlotsPerPage)))
                } else {
                    vecidxPlots <- seq((p * scaNPlotsPerPage + 1), scaNIDs)
                }
                print(plot_grid(plotlist = lsgplotsID[vecidxPlots], align = "h", 
                                nrow = scaNPlotsPerPage/2, ncol = 2, 
                                rel_widths = c(1, 1), rel_heights = c(1, 1, 1)))
            }
        } else {
            pdf(dirFileOut)
            for (p in seq(1, scaNIDs)) {
                print(lsgplotsID[[p]])
            }
        }
        dev.off()
        graphics.off()
    }
    
    # Return raw list of gplots
    return(lsgplotsID)
}

### Return model fits

#' Return a gene x sample and gene x time point (for each condition) matrix with model fits.
#' 
#' @seealso Called seperately by user.
#' 
#' @param objectImpulseDE2 (instance of class ImpulseDE2Object)
#' ImpulseDE2 output object to create heatmap from.
#' 
#' @return 
#' \itemize{
#' \item sample (data.frame)
#' Data frame with model fit by gene and sample.
#' \item case (data.frame)
#' Data frame with model fit by gene and time point in case condition
#' \item control (data.frame)
#' Data frame with model fit by gene and time point in control condition
#' }
#' 
#' @examples
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 0,
#' scaNImp          = 50,
#' scaNLin          = 0,
#' scaNSig          = 50)
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' boolBeta2       = FALSE,
#' vecCovFactor    = NULL,
#' vecCovContinous = NULL,
#' boolIdentifyTransients = FALSE,
#' scaNProc        = 1 )
#' modelFits <- computeModelFits(objectImpulseDE2=objectImpulseDE2)
#' head(modelFits$case)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
computeModelFits <- function(objectImpulseDE2){
    
    dfAnnot <- get_dfDESeqAnnotationProc(obj=objectImpulseDE2)
    
    # Compute by time point and by sample:
    vecTimePointsToEval <- dfAnnot$Time
    # Get sample group assignment indices: Time and batch
    vecTimepointsUnique <- sort(unique(vecTimePointsToEval))
    vecidxTimepoint <- match(vecTimePointsToEval, vecTimepointsUnique)
    # Get size factors for all samples
    vecSfValues <- get_vecSizeFactors(obj=objectImpulseDE2)
    
    calcImpulseValueCondition <- function(obj, condition, dfAnnot){
        # Gene x sample matrix for samples in "condition"

        dfAnnotation <- dfAnnot[dfAnnot == condition]
        samples <- dfAnnotation$Sample

        # Time points:
        vecTimepoints <- dfAnnotation$Time
        vecTimepointsUnique <- sort(unique(vecTimepoints))
        vecidxTimepoint <- match(vecTimepoints, vecTimepointsUnique)

        # Scaling factors 
        # Size factors (per-sample):
        vecSfValues <- get_vecSizeFactors(obj=objectImpulseDE2)[samples]
        # Correction factors (per-condition):
        # Names are covariates 
        vecCorrectionFactors <- get_lsModelFits(obj=objectImpulseDE2)[[condition]][[x]]$lsImpulseFit$vecCorrectionFactors

        # Now do things on a gene-level:
        matImpulseValueSample <- do.call(rbind, lapply(get_vecAllIDs(obj=objectImpulseDE2), function(x) {
            # Get impulse values for samples in "condition"
            # One value per sample 
            vecImpulseValues <- ImpulseDE2:::evalImpulse_comp(
                vecImpulseParam = get_lsModelFits(obj=objectImpulseDE2)[[condition]][[x]]$lsImpulseFit$vecImpulseParam, 
                vecTimepoints = vecTimepointsUnique)[vecidxTimepoint] 
            # Make sure samples match up 
            stopifnot(names(vecImpulseValues) == names(vecSfValues))
            # Multiply impulse values by size factors
            vecModelValues <- vecImpulseValues * vecSfValues    
            # Multiply model values by correction factors from covariates
            if(!is.null(vecCorrectionFactors)){
                vecModelValues <- vecModelValues * vecCorrectionFactors
            }
            return(vecModelValues)
            })
        )

        rownames(matImpulseValueSample) <- get_vecAllIDs(obj=objectImpulseDE2)
        colnames(matImpulseValueSample) <- samples
        return(matImpulseValueSample)
    }

    # Get impulse values for "case" samples
    matImpulseValueSample <- calcImpulseValueCondition(obj=objectImpulseDE2,
        condition="case", dfAnnot=dfAnnot)
    if(get_boolCaseCtrl(obj=objectImpulseDE2)){
        # Get impulse values for "control" samples
        matImpulseValueSampleControl <- calcImpulseValueCondition(obj=objectImpulseDE2,
            condition="control", dfAnnot=dfAnnot)
        # Concatenate matrices
        matImpulseValueSample <- cbind(matImpulseValueSample, matImpulseValueSampleControl)
    }

    # matImpulseValueSample <- do.call(rbind, lapply(get_vecAllIDs(obj=objectImpulseDE2), function(x) {
    #     # Get case impulse values for gene x
    #     vecImpulseValuesCase <- ImpulseDE2:::evalImpulse_comp(
    #         vecImpulseParam = get_lsModelFits(obj=objectImpulseDE2)$case[[x]]$lsImpulseFit$vecImpulseParam, 
    #         vecTimepoints = vecTimepointsUnique)[vecidxTimepoint]
    #     vecModelValues <- vecImpulseValuesCase * vecSfValues
    #     if (!is.null(objectImpulseDE2@lsModelFits$IdxGroups$case$lsvecidxBatch)) {
    #         vecBatchFactorsCase <- get_lsModelFits(
    #             obj=objectImpulseDE2)$case[[x]]$lsImpulseFit$lsvecBatchFactors
    #         for (i in seq(length(vecBatchFactorsCase))) {
    #             vecModelValues <- vecModelValues * vecBatchFactorsCase[[i]][
    #                 objectImpulseDE2@lsModelFits$IdxGroups$case$lsvecidxBatch[[i]]
    #             ]
    #         }
    #     }
        
    #     if (!is.null(objectImpulseDE2@lsModelFits$IdxGroups$control)) {
    #         vecImpulseValuesCtrl <- ImpulseDE2:::evalImpulse_comp(
    #             vecImpulseParam = 
    #                 get_lsModelFits(obj=objectImpulseDE2)$control[[x]]$
    #                 lsImpulseFit$vecImpulseParam, 
    #             vecTimepoints = dfAnnot$Time
    #         )
    #         vecModelValuesCtrl <- vecImpulseValuesCtrl * vecSfValues
    #         if (!is.null(objectImpulseDE2@lsModelFits$IdxGroups$case$lsvecidxBatch)) {
    #             vecBatchFactorsCtrl <- get_lsModelFits(
    #                 obj=objectImpulseDE2)$control[[x]]$lsImpulseFit$lsvecBatchFactors
    #             for (i in seq(length(vecBatchFactorsCase))) {
    #                 vecModelValuesCtrl <- vecModelValuesCtrl * vecBatchFactorsCtrl[[i]][
    #                     objectImpulseDE2@lsModelFits$IdxGroups$control$lsvecidxBatch[[i]]
    #                 ]
    #             }
    #         }
    #         vecModelValues[dfAnnot$Condition == "control"] <- 
    #             vecModelValuesCtrl[dfAnnot$Condition == "control"]
    #     }
    #     return(vecModelValues)
    # }))
    # rownames(matImpulseValueSample) <- dfAnnot$vecAllIDs
    # colnames(matImpulseValueSample) <- dfAnnot$Sample

    # Compute by time point in case condition:
    matImpulseValueCase <- do.call(rbind, lapply(objectImpulseDE2@vecAllIDs, function(x) {
        evalImpulse_comp(
            vecImpulseParam = get_lsModelFits(obj=objectImpulseDE2)$case[[x]]$lsImpulseFit$vecImpulseParam, 
            vecTimepoints = sort(unique(dfAnnot$Time), decreasing = FALSE)
        )
    }))
    rownames(matImpulseValueCase) <- objectImpulseDE2@vecAllIDs
    colnames(matImpulseValueCase) <- sort(unique(dfAnnot$Time), decreasing = FALSE)
    
    # Compute by time point in control condition:
    if (objectImpulseDE2@boolCaseCtrl) {
        matImpulseValueCtrl <- do.call(rbind, lapply(objectImpulseDE2@vecAllIDs, function(x) {
            evalImpulse_comp(
                vecImpulseParam = get_lsModelFits(obj=objectImpulseDE2)$control[[x]]$lsImpulseFit$vecImpulseParam, 
                vecTimepoints = sort(unique(dfAnnot$Time), decreasing = FALSE)
            )
        }))
        rownames(matImpulseValueCtrl) <- objectImpulseDE2@vecAllIDs
        colnames(matImpulseValueCtrl) <- sort(unique(dfAnnot$Time), decreasing = FALSE)
    } else {
        matImpulseValueCtrl <- NULL
    }
    
    return(list(sample = matImpulseValueSample, 
                case = matImpulseValueCase, 
                control = matImpulseValueCtrl))
}

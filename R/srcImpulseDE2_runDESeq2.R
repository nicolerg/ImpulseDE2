### Estimate dispersions with DESeq2

#' Wrapper function for running DESeq2
#' 
#' Run DESeq2 and extract dispersion parameter estimates.
#' Catch and remove dispersion outlier exception on samples
#' with zero-count observations.
#' 
#' @seealso Called by \link{runImpulseDE2}.
#' 
#' @param dfAnnotationProc (data frame samples x covariates) 
#' {Sample, Condition, Time (numeric), TimeCateg (str)
#' (and covariates if given).}
#' Processed annotation table with covariates for each sample.
#' @param matCountDataProc (matrix genes x samples)
#' Read count data.
#' @param boolCaseCtrl (bool) 
#' Whether to perform case-control analysis. Does case-only
#' analysis if FALSE.
#' @param vecCovFactor (vector of strings number of categorical covariates)
#' Categorial covariates to adjust for.
#' Names refer to columns in dfAnnotation.
#' @param vecCovContinuous (vector of strings number of continuous covariates)
#' Continuous covariates to adjust for.
#' Names refer to columns in dfAnnotation. 
#'		
#' @return (numeric vector length number of genes)
#' Dispersion parameter estimates for each gene.
#' In format of parameter size of \link{dnbinom}
#' which is 1/dispersion factor of DESeq2.
#'    
#' @import DESeq2
#' @importFrom S4Vectors mcols
#' 
#' @author David Sebastian Fischer
runDESeq2 <- function(
    dfAnnotationProc, 
    matCountDataProc,
    boolCaseCtrl,
    vecCovFactor,
    vecCovContinuous){

    # Checks whether a model matrix from a data.frame is full rank
    isFullRank <- function(df){
        matModelMatrix = model.matrix(~., data=df)
        return(ncol(matModelMatrix) == rankMatrix(matModelMatrix))
    }

    vecCovAll = c(vecCovContinuous, vecCovFactor)
    
    # Get gene-wise dispersion estimates
    # var = mean + alpha * mean^2, alpha is dispersion
    # DESeq2 dispersion is 1/size used dnbinom (used in cost function
    # for evaluation of likelihood)
    
    dds <- NULL
    if(!boolCaseCtrl){
        # Case-only
        if(is.null(vecCovAll)){
            # No adjustment
            dds <- suppressWarnings( DESeqDataSetFromMatrix(
                countData = matCountDataProc,
                colData = dfAnnotationProc,
                design = ~ TimeCateg) )
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersions(dds)
            
        } else {
            # Adjustment 
            cont = paste0(c("~ TimeCateg", vecCovAll), 
                collapse = " + ")
            tryCatch({
                print(cont)
                dds <- suppressWarnings( DESeqDataSetFromMatrix(
                    countData = matCountDataProc,
                    colData = dfAnnotationProc,
                    design = eval(parse(text=cont)) ) )
                dds <- estimateSizeFactors(dds)
                dds <- estimateDispersions(dds)
            }, error=function(strErrorMsg){
                print(strErrorMsg)
                contred = paste0("~ ", paste0(vecCovAll, collapse = " + "))
                print(paste0(
                    "WARNING: DESeq2 failed on full model ",
                    "- dispersions may be inaccurate.",
                    "Estimating dispersions on reduced model ",
                    "formulation [full = ", contred,
                    ", reduced = ~1]. Supply externally generated",
                    " dispersion parameters via ",
                    "vecDispersionsExternal if there is a more ",
                    "accurate model for your data set."))
                warning("Warning generated in dispersion factor",
                        " estimation, read stdout.")
            }, finally={
                if(is.null(dds)){
                    contred = paste0("~ ", paste0(vecCovAll, collapse = " + "))
                    dds <- suppressWarnings( DESeqDataSetFromMatrix(
                        countData = matCountDataProc,
                        colData = dfAnnotationProc,
                        design = eval(parse(text=contred)) ) )
                    dds <- estimateSizeFactors(dds)
                    dds <- estimateDispersions(dds)
                }
            })
        } 
    } else {
        # Case-control
        if(is.null(vecCovAll)){
            # No correction
            # Catch non full-rank design matrix 
            # e.g. conditions have mutually 
            #exclusive sets of timepoints
            tryCatch({
                dds <- suppressWarnings( DESeqDataSetFromMatrix(
                    countData = matCountDataProc,
                    colData = dfAnnotationProc,
                    design = ~Condition + Condition:TimeCateg) )
                dds <- estimateSizeFactors(dds)
                dds <- estimateDispersions(dds)
            }, error=function(strErrorMsg){
                print(strErrorMsg)
                print(paste0(
                    "WARNING: DESeq2 failed on full model ",
                    "- dispersions may be inaccurate.",
                    "Estimating dispersions on reduced model",
                    " formulation [full = ~Condition,",
                    " reduced = ~1]. Supply externally ",
                    "generated dispersion parameters via ",
                    "vecDispersionsExternal if there is a more",
                    " accurate model for your data set."))
                warning("Warning generated in dispersion factor",
                        " estimation, read stdout.")
            }, finally={
                if(is.null(dds)){
                    dds <- suppressWarnings( DESeqDataSetFromMatrix(
                        countData = matCountDataProc,
                        colData = dfAnnotationProc,
                        design = ~Condition) )
                    dds <- estimateSizeFactors(dds)
                    dds <- estimateDispersions(dds)
                }
            })
            
        } else {
            # With correction 
            # Catch non full-rank design matrix 
            # This is also checked for during the preprocessing step 
            cont = paste0(c("~ Condition + Condition:TimeCateg", vecCovAll), collapse = " + ")
            tryCatch({
                dds <- suppressWarnings( DESeqDataSetFromMatrix(
                    countData = matCountDataProc,
                    colData = dfAnnotationProc,
                    design = eval(parse(text=cont)) ) )
                dds <- estimateSizeFactors(dds)
                dds <- estimateDispersions(dds)
            }, error=function(strErrorMsg){
                print(strErrorMsg)
                print(paste0(
                    "WARNING: DESeq2 failed on full model ",
                    "- dispersions may be inaccurate.",
                    "Estimating dispersions on reduced model.",
                    " Supply externally generated dispersion parameters via ",
                    "vecDispersionsExternal if there is a more accurate ",
                    "model for your data set."))
                warning("Warning generated in dispersion factor estimation.",
                        " Read stdout.")
            }, finally={
                if(is.null(dds)){
                    # One of three cases:
                    # 1. Covariates are linear combinations of time points (covariates + time is not full rank)
                    # 2. Conditions (case and control) have mutually exclusive sets of time points (condition + time is not full rank)
                    # 3. Covariates are linear combinations of other covariates (covariates alone are not full rank)
                    # 4. Covariates are linear combinations of conditions (covariates + condition is not full rank)
                    if(!isFullRank(dfAnnotationProc[,c("TimeCateg", vecCovAll)])){
                        # Case 1
                        contred = paste0(c("~ Condition", vecCovAll), collapse = " + ")
                        print(paste0("Model matrix based on covariates {",
                               paste(vecCovAll, collapse=','), 
                               "} and Time is not full rank: ",
                               "There are covariates which are linear",
                               " combinations of the time points."))
                        print(paste0("Using reduced model formulation ",
                               "[full= ", contred, ", ",
                               "reduced= ~1]."))
                        dds <- suppressWarnings( DESeqDataSetFromMatrix(
                            countData = matCountDataProc,
                            colData = dfAnnotationProc,
                            design = eval(parse(text=contred)) ) )
                        dds <- estimateSizeFactors(dds)
                        dds <- estimateDispersions(dds)
                    }else if(!isFullRank(dfAnnotationProc[,c("Condition", "TimeCateg")])){
                        # Case 2
                        contred = paste0(c("~ Condition", vecCovAll), collapse = " + ")
                        print(paste0("Model matrix based on Condition",
                               " and Time is not full rank: ",
                               "Conditions case and control have",
                               " mutually exclusive sets of timepoints."))
                        print(paste0("Using reduced model formulation ",
                               "[full= ", contred, ", ",
                               "reduced= ~1]."))
                        dds <- suppressWarnings( DESeqDataSetFromMatrix(
                            countData = matCountDataProc,
                            colData = dfAnnotationProc,
                            design = eval(parse(text=contred)) ) )
                        dds <- estimateSizeFactors(dds)
                        dds <- estimateDispersions(dds)
                    }else if(!isFullRank(dfAnnotationProc[,vecCovAll])){
                        # Case 3
                        contred = "~ Condition + Condition:TimeCateg"
                        print(paste0("Model matrix based on",
                               " covariates {",paste(vecCovAll, collapse=','), "} is not full rank: ",
                               "Some covariates are linear combinations of each other."))
                        print(paste0("Using reduced model formulation ",
                               "[full= ", contred, ", ",
                               "reduced= ~1]."))
                        dds <- suppressWarnings( DESeqDataSetFromMatrix(
                            countData = matCountDataProc,
                            colData = dfAnnotationProc,
                            design = eval(parse(text=contred)) ) )
                        dds <- estimateSizeFactors(dds)
                        dds <- estimateDispersions(dds)
                    }else if(!isFullRank(dfAnnotationProc[,c("Condition", vecCovAll)])){
                        # Case 4
                        contred = paste0(c("~ TimeCateg", vecCovAll), collapse = " + ")
                        print(paste0("Model matrix based on covariates {",
                               paste(vecCovAll, collapse=','), 
                               "} and Condition is not full rank: ",
                               "There are covariates which are linear",
                               " combinations of Condition (case/control)."))
                        print(paste0("Using reduced model formulation ",
                               "[full= ", contred, ", ",
                               "reduced= ~1]."))
                        dds <- suppressWarnings( DESeqDataSetFromMatrix(
                            countData = matCountDataProc,
                            colData = dfAnnotationProc,
                            design = eval(parse(text=contred)) ) )
                        dds <- estimateSizeFactors(dds)
                        dds <- estimateDispersions(dds)
                    }else{
                        stop(paste("Congratulations, I don't know how you got here.", 
                            "Please check previous warning and error messages.",
                            "There is unexpected confounding between your covariates {",
                            paste(vecCovAll, collapse=','),"}, Time,",
                            "and Condition that we did not account for."))
                    }
                }
            })
        }
    }
    
    vecDispersionsInv <- mcols(dds)$dispersion
    # Catch dispersion trend outliers at the upper boundary
    # (alpha = 20 ->large variance)
    # which contain zero measurements: 
    # The zeros throw off the dispersion estimation 
    # in DESeq2 which may converge to the upper bound 
    # even though the obesrved variance is small.
    # Avoid outlier handling and replace estimates by 
    # MAP which is more stable in these cases.
    # Note: the upper bound 20 is not exactly reached - 
    # there is numeric uncertainty here - use > rather than ==
    vecindDESeq2HighOutliesFailure <- !is.na(mcols(dds)$dispOutlier) & 
        mcols(dds)$dispOutlier==TRUE &
        mcols(dds)$dispersion>(20-10^(-5)) & 
        apply(matCountDataProc, 1, function(gene) any(gene==0) )
    vecDispersionsInv[vecindDESeq2HighOutliesFailure] <- 
        mcols(dds)$dispMAP[vecindDESeq2HighOutliesFailure]
    if(sum(vecindDESeq2HighOutliesFailure)>0){
        print(paste0("Corrected ", sum(vecindDESeq2HighOutliesFailure),
                     " DESEq2 dispersion estimates which ",
                     "to avoid variance overestimation and loss of ",
                     "discriminatory power for model selection."))
    }
    # DESeq2 uses alpha=1/phi as dispersion
    vecDispersions <- 1/vecDispersionsInv
    names(vecDispersions) <- rownames(dds)
    
    return(vecDispersions)
}